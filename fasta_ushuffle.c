/*
   fasta_ushuffle - Shuffles sequences in a FASTA file.
   Copyright (C) 2010 Assaf Gordon (gordon@cshl.edu)

   Released under the same license as uShuffle (see below):

=== uShuffle Copyright notice ===
 * Copyright (c) 2007
 *   Minghui Jiang, James Anderson, Joel Gillespie, and Martin Mayne.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice,
 *      this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *      this list of conditions and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 * 3. The names of its contributors may not be used to endorse or promote
 *      products derived from this software without specific prior written
 *      permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 *	fasta_ushuffle.c - command-line interface of uShuffle
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <err.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <stdbool.h>
#include "ushuffle.h"

//Hard-coded limits, seems resonable for next-gen (short) reads
#define	MAX_ID_SIZE 32678
#define MAX_SEQUENCE_SIZE (1000000)

#define VERSION "0.2"

#define HELPTEXT \
"fasta_ushuffle: shuffles biological sequences while preserving the k-let counts.\n" \
"\n" \
"VERSION " VERSION "\n" \
"\n" \
"\nCopyright (C) 2010 A. gordon (gordon@cshl.edu).\n" \
"\nUses the uShuffle library code by: Minghui Jiang, James Anderson, Joel Gillespie, and Martin Mayne.\n" "\n" \
"Usage: fasta_ushuffle [-r N] [-h] [-o] [-n N] [-k N] [-s N] < INPUT.FA > OUTPUT.FA\n" \
"\n" \
" -h		This help screen\n" \
" -o            Print original (unshuffled) in output file.\n" \
" -k N		specifies the let size\n" \
" -s N		specifies the seed for random number generator.\n" \
" -n N          For each input sequence, print N permutations (default is 1).\n" \
"               Use this only for debugging.\n" \
" -r N          Retry N times to find a new shuffle (Default is 10). After N retries, a warning is printed, and a non-shuffled sequence will be written.\n" \
"\n" \
"Nucleotide sequences in the input FASTA file must be in a single line.\n" \
"This is a valid input file:\n" \
"  >dummy1\n" \
"  AGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGAGTG\n" \
"  >dummy2\n" \
"  CTGAGAGTCACACATGATTTTACAACAACCATGAAG\n" \
"\n" \
"This is not a valid input file:\n" \
"  >dummy1\n" \
"  AGTAGTAGTAGTAGTAGTAGTAGTAG\n" \
"  TAGTAGAGTG\n" \
"  >dummy2\n" \
"  CTGAGAGTCACACATGATTTTACAAC\n" \
"  AACCATGAAG\n" \
"\n" \
"Use fasta_formatter (from the FASTX-Toolkit) to re-format a multiline fasta file.\n" \
"\n"

void showhelp()
{
	fprintf(stderr, HELPTEXT );
	exit(0);
}

bool is_valid_nucleotide_string(const char *s)
{
	if (s==NULL)
		return false;
	if (s[0]==0)
		return false;

	while ( *s != 0 ) {
		switch(*s)
		{
		case 'A':
		case 'C':
		case 'G':
		case 'T':
		case 'R':
		case 'Y':
		case 'S':
		case 'W':
		case 'K':
		case 'M':
		case 'B':
		case 'D':
		case 'H':
		case 'V':
		case 'N':
		case 'a':
		case 'c':
		case 'g':
		case 't':
		case 'r':
		case 'y':
		case 's':
		case 'w':
		case 'k':
		case 'm':
		case 'b':
		case 'd':
		case 'h':
		case 'v':
		case 'n':
			break;

		default:
			return false;
		}

		++s;
	}
	return true;
}


/*
   Poor man's FASTA parser and validator.

   Reads two lines from STDIN, validates them as FASTA format.
 */
bool read_fasta_record(char* /*OUTPUT*/ fasta_id, int max_id_size,
			char * /*output*/ fasta_sequence, int max_sequence_size,
			unsigned long line)
{
	if (fgets(fasta_id,max_id_size,stdin)==NULL)
		return false; //EOF - this is not an error

	//
	// First line - FASTA ID
	//

	//The input line was too long, and was trimmed -  Abort the program.
	//We don't do dynamic memory allocation here for simplicty.
	size_t id_len = strlen(fasta_id);
	if (id_len==max_id_size-1) {
		fprintf(stderr,"Internal error: got a too-long input line (line %lu). Please incease the value of MAX_ID_SIZE (currently = %d)\n", line, MAX_ID_SIZE);
		exit(1);
	}

	//Too short ? not a valid FASTA identifier
	if (id_len<2) {
		fprintf(stderr,"Input error: got too-short ID line (line %lu).\n", line);
		exit(1);
	}

	//Chomp the ID line
	if (fasta_id[id_len-1]=='\n')
		fasta_id[id_len-1]=0;

	//FASTA identifiers must begin with '>'
	if (fasta_id[0]!='>') {
		if (is_valid_nucleotide_string(fasta_id)) {
			//A Multiline FASTA file - detect and warn the user
			fprintf(stderr,"Input error: input looks like a multi-line FASTA file (line %lu should start with '>' but contains nucleotide sequence). This program requires a single-line FASTA file. Use 'fasta_formatter' to re-format the input file.\n", line);
			exit(1);
		}

		//Otherwise - just complain
		fprintf(stderr,"Input error: Invalid FASTA identifier on line %lu (expecting line with '>').\n", line);
		exit(1);
	}

	/********************************
	 * the nuceleotide sequence line
	 *************************************/
	++line;
	if (fgets(fasta_sequence,max_sequence_size,stdin)==NULL) {
		fprintf(stderr,"Error: Missing nucleotide sequence line in input FASTA file (line %lu\n", line);
		exit(1);
	}

	//The input line (sequence line) was too long, and was trimmed -  Abort the program.
	//We don't do dynamic memory allocation here for simplicty.
	if (strlen(fasta_sequence)==max_sequence_size-1) {
		fprintf(stderr,"Internal error: got a too-long input line (line %lu). Please incease the value of MAX_SEQUENCE_SIZE (currently = %d)\n", line, max_sequence_size);
		exit(1);
	}

	//Valid nucleotide string?
	//Chomp the ID line
	size_t seq_len = strlen(fasta_sequence);

	//chomp
	if (fasta_sequence[seq_len-1]=='\n')
		fasta_sequence[seq_len-1]=0;

	if (!is_valid_nucleotide_string(fasta_sequence)) {
		fprintf(stderr,"Input error: Invalid input file, expecting nucleotide sequence line on line %lu\n", line);
		exit(1);
	}

	return true;
}

void print_shuffle_sequence_perm(int k, int permutations_count, const char*id, const char*sequence)
{
	int l;
	char *t=NULL;
	int i;

	l = strlen(sequence);
	if ((t = malloc(l + 1)) == NULL)
		err(1,"malloc failed");
	t[l] = '\0';

	shuffle1(sequence, l, k);
	for (i = 0; i < permutations_count; i++) {
		shuffle2(t);
		printf("%s-perm%d\n", id, i+1);
		printf("%s\n", t);
	}
	shuffle_reset();

	free(t);
}

void print_shuffle_sequence_retries(int k, int retries_count, const char*id, const char*sequence)
{
	int l;
	char *t=NULL;
	int i;

	l = strlen(sequence);
	if ((t = malloc(l + 1)) == NULL)
		err(1,"malloc failed");
	t[l] = '\0';

	shuffle1(sequence, l, k);

	i = 0 ;
	while ( i < retries_count ) {
		shuffle2(t);
		if (strncmp(sequence, t, l) != 0) {
			printf("%s\n", id);
			printf("%s\n", t);
			break;
		}
		i++;
	}
	if (i>=retries_count) {
		fprintf(stderr,"WARNING: failed to find new shuffle for sequence \"%s\" (%s) after %d retries\n", id, sequence, retries_count);
		printf("%s\n", id);
		printf("%s\n", t);
	}
	shuffle_reset();

	free(t);
}

int main(int argc, char **argv)
{
	char *s = NULL, *t;
	int n = 1, k = 2, b = 0;
	struct rusage r1, r2;
	struct timeval t1, t2;
	struct timeval tv;
	double u1;
	unsigned long seed;
	int i, l;
	int c;
	unsigned long line=1;
	bool show_original=false;
	int max_retries=10;

	char*	fasta_id;
	char*	fasta_sequence;

	if ((fasta_id = malloc(MAX_ID_SIZE))==NULL)
		err(1,"malloc(%d) failed", MAX_ID_SIZE);
	if ((fasta_sequence = malloc(MAX_SEQUENCE_SIZE))==NULL)
		err(1,"malloc(%d) failed", MAX_SEQUENCE_SIZE);

	gettimeofday(&tv, NULL);
	seed = (unsigned long) tv.tv_sec;

	// Parse command line options
	while ( (c=getopt(argc, argv, "ok:n:s:hr:"))!=-1) {
		switch (c)
		{
		case 'o':
			show_original = true;
			break;

		case 'n':
			n = atoi(optarg);
			if (n<=0) {
				fprintf(stderr,"Error: invalid -n value (%s). Must be a number larger than zero.", optarg);
				exit(1);
			}
			break;
		case 'k':
			k = atoi(optarg);
			if (k<=0) {
				fprintf(stderr,"Error: invalid -k value (%s). Must be a number larger than zero.", optarg);
				exit(1);
			}
			break;

		case 'r':
			max_retries = atoi(optarg);
			if (max_retries<=0) {
				fprintf(stderr,"Error: invalid -r value (%s). Must be a number larger than zero.", optarg);
				exit(1);
			}
			break;

		case 's':
			seed = atoi(optarg);
			break;

		default:
		case 'h':
			showhelp();
		}
	}

	srandom(seed);
	set_randfunc((randfunc_t) random);

	while (read_fasta_record(fasta_id,MAX_ID_SIZE, fasta_sequence, MAX_SEQUENCE_SIZE, line)) {
		line+=2;

		if (show_original) {
			printf("%s-unshuffled\n", fasta_id);
			printf("%s\n", fasta_sequence);
		}

		if (n>1) {
			print_shuffle_sequence_perm(k, n, fasta_id, fasta_sequence);
		} else {
			print_shuffle_sequence_retries(k, max_retries, fasta_id, fasta_sequence);
		}
	}

	free(fasta_id);
	free(fasta_sequence);

	return 0;
}
