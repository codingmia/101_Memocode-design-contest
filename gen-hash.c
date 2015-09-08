/** 
	gen-hash.c
	 Creates a dynamic hash table structure, based on 16-bit 
	 binary strings from the supplied genome_file. 
	 Converts the dynamic hash structure into two static structures, and
	 outputs those as two binary files. 
	 Prefix table and Hash table
 **/

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#include <string.h>
#include <stdint.h>

/*
	#define HASH_LENGTH 2 // 2 Bytes = 16 bits
	#define TABLE_LENGTH 65536 // 2^ 16
	hash_entry_t pointer: a list of hash_entry
	hash_entry_t pointer pointer: a list of list of hash entry
*/

#define HASH_LENGTH 3
#define TABLE_LENGTH 16777216 //24bits: there are 2^24 variations of pairs

#define NUM_BLOCKS 8

/* Dynamic hash table entry */
typedef struct hash_entry_s {
  uint32_t position;
  struct hash_entry_s *next;	
} hash_entry_t; //alias

hash_entry_t **hash_heads, **hash_tails;

void gen_hash(uint8_t *reference_base, off_t reference_total, uint32_t seq_min, uint32_t seq_max) {
	
  uint32_t i, j, k, index;
  uint32_t seq, tseq;
  uint8_t bytes[25], skip;
  hash_entry_t *hash_tmp;
	
  /* Allocate memory for the hash table */
  hash_heads = (hash_entry_t **)calloc(TABLE_LENGTH, sizeof(hash_entry_t *));
  hash_tails = (hash_entry_t **)calloc(TABLE_LENGTH, sizeof(hash_entry_t *));


  /* As we loop through the reference genome, we can stop indexing when
   * there is no more room for a full match anyway 
   * reference_total - 25 to reference_total still has one full match
   */
  for (i = 0; i < reference_total-24; i++) {
    /* bytes[25] stores one 100 bp 
     *(AGCT, need to 2 bits to encode one base, so 100 bases will need 200bits  = 25 Bytes)
     *Calculate 3 byte length prefix
     */
    for (k = 0; k < HASH_LENGTH+1; k++) {
      bytes[k] = reference_base[i+k];
    }
    //Processed - One byte, four bases - AGCT
    //Processed - sequence start with T ---, C ---, G ---, A ---
    for (j = 0; j < 4; j++) {
      /* 
      	The prefix tree is a function of HASH_LENGTH. Look at that many bytes
      	seq stores the prefix by first 12 bp of 100 bp, namely, 3 Bytes of 25 Bytes
       */

      seq = 0; 
      for (k = 0; k < HASH_LENGTH; k++) {
		seq |= bytes[k] << ((HASH_LENGTH-k-1)*8);//16,  8 , 0
      } // x << 16 + x << 8 + x is 24 bits

      /* Try to find some sequences we can skip. In this version, we skip
       *  every sequence that is all 0s*/

      skip = 0;
      //first three bytes prefix are zeros, continue checking the remaining.
      if (seq == 0) {
		tseq = 0;
		for (k = HASH_LENGTH+1; k < 24; k++) {
		  tseq |= reference_base[i+k];
		}
		if (tseq == 0) {
		  //fprintf(stderr, "Skipped a reference sequence\n");
		  skip = 1;
		}
      }
      
      /* Partition the hash table based on seq_min and seq_max */
      if (seq >= seq_min && seq < seq_max && skip == 0) {

		//printf("Calculated sequence @ %i of %04x\n", 4*i+j,seq);
	    
		/* Allocate a new hash node */
		hash_tmp = (hash_entry_t *)malloc(sizeof(hash_entry_t));
		if (!hash_tmp) {
		  perror("Error when allocating hash table memory\n");
		}

		hash_tmp->position = 4*i+j;//actual location of the current prefix in the bytes array.
		hash_tmp->next = NULL;
	     
	    //Place the current prefix in the Prefix table
		/* If it's the first node in the row, set the pointers appropriately */
		if (hash_heads[seq] == NULL) {
		  //printf("Found first entry for %06x, at position %i\n", seq, 4*i+j);
		  hash_heads[seq] = hash_tmp; 
		}
		else {
		  hash_tails[seq]->next = hash_tmp;
		}
		hash_tails[seq] = hash_tmp;//append hash_tmp to the tail of hash entry list
      }

      /* Update the sequence value by shifting out 2 bits, and shifting
       * in 2 bits from the next reference_base byte. The fasta/fastq 
       * formats are read from right to left for each byte. */
      update_seq:
      for (k = 0; k < HASH_LENGTH+1; k++) {
		bytes[k] = bytes[k] >> 2;
		bytes[k] = bytes[k] & 0x3f;
		bytes[k] = bytes[k] | ((bytes[k+1] & 0x3) << 6);
      }
    }//end of j [0, 3]
  }
  return;	
}

void bin_hash(uint32_t seq_min, uint32_t seq_max) {
  FILE *fp1, *fp2;
  uint32_t i, j, matches, prev_matches;
  hash_entry_t *hash_ptr, *prev_ptr;
  uint32_t *table1, *table2;
  off_t table2_total;
  struct stat file_status;
  
  if (seq_min == 0) {
    fp1 = fopen("hash_table1.bin", "wb");
  }
  else {
    fp1 = fopen("hash_table1.bin", "ab");
  }

  if (fp1 == NULL) {
    fprintf(stderr, "Error: could not open hash_table1.bin for writing\n");
    return;		
  }


  if (seq_min == 0) {
    fp2 = fopen("hash_table2.bin", "wb");
  }
  else {
    fp2 = fopen("hash_table2.bin", "ab");
  }  
  if (fp2 == NULL) {
    fprintf(stderr, "Error: could not open hash_table2.bin for writing\n");
    return;		
  }

  if (fstat(fileno(fp2), &file_status)) {
    fprintf(stderr, "Error checking hash_table2.bin\n");
    perror((const char *) 0);
  }
  table2_total = file_status.st_size;


  table1 = (uint32_t *)malloc((seq_max-seq_min)*sizeof(uint32_t));
  if (!table1) {
    perror("Error when allocating hash table memory\n");
  }

  /* The start point of the current hash_table1 needs to be the size of 
   * the current hash_table2 */
  table1[0] = table2_total/4;
  prev_matches = 0;
  for (i = seq_min; i < seq_max; i++) {

    /* First, find the length of each row, and dump to the file */
    matches = 0;
    for (hash_ptr = hash_heads[i]; hash_ptr != NULL; 
	 hash_ptr = hash_ptr->next) {
      matches++;
    }

    //printf("Found %d matches for seq %06x\n", matches, i);
    table2 = (uint32_t *)malloc(matches*sizeof(uint32_t));
    if (!table2) {
      perror("Error when allocating hash table memory\n");
    }
    matches = 0;
    prev_ptr = NULL;
    for (hash_ptr = hash_heads[i]; hash_ptr != NULL; 
	 hash_ptr = hash_ptr->next) {
      table2[matches++] = hash_ptr->position;

      if (prev_ptr != NULL) {
	free(prev_ptr);
      }
      prev_ptr = hash_ptr;
    }


    fwrite(table2, sizeof(uint32_t), matches, fp2);
    free(table2);
    
    /* Table1 is easy, first entry is 0 and then we increment by matches */
    if (i != seq_min) {
      table1[i-seq_min] = table1[i-seq_min-1] + prev_matches;
    }
    prev_matches = matches;

  }
  
  fwrite(table1, sizeof(uint32_t), seq_max-seq_min, fp1);

  free(table1);
  fclose(fp1);
  fclose(fp2);


  /* Free up some memory */
  free(hash_heads);
  free(hash_tails);

  
  return;
}

/*
 * Process command-line arguments, map reference and sequence data into
 * memory, and call hash()
 */
int main(int argc, const char *argv[]) {

  const uint8_t *reference_filename;
  int32_t reference_fd = -1;
  off_t reference_total;
  int32_t i;

  long page_size;
  void *reference_base;
  struct stat file_status;  

  page_size = sysconf(_SC_PAGE_SIZE); /* needed for mmap */

  if (argc < 2) goto usage;

  reference_filename = argv[1];

  if ((reference_fd = open(reference_filename, O_RDONLY)) < 0) {
    fprintf(stderr, "Error opening reference file \"%s\": ",
	    reference_filename);
    perror((const char *) 0);
    goto usage;
  }

  if (fstat(reference_fd, &file_status)) {
    fprintf(stderr, "Error checking reference file \"%s\": ",
	    reference_filename);
    perror((const char *) 0);
    goto usage;
  }

  reference_total = file_status.st_size;

  /* mmap the reference data */
  reference_base = mmap( (void *) 0, reference_total, PROT_READ, MAP_SHARED, 
  		reference_fd, 0);

  if (reference_base == MAP_FAILED) {
    perror("Error when attempting to map the reference file");
    goto closeexit;
  }

  /*NUM_BLOCKS = 8
  	TABLE_LENGTH/NUM_BLOCKS = 2^24 / 8
  	*/
  for (i = 0; i < NUM_BLOCKS; i++) {
    gen_hash(reference_base, reference_total, TABLE_LENGTH/NUM_BLOCKS*i, TABLE_LENGTH/NUM_BLOCKS*(i+1));
    bin_hash(TABLE_LENGTH/NUM_BLOCKS*i, TABLE_LENGTH/NUM_BLOCKS*(i+1));
  }

 unmap_references:
  if (munmap(reference_base, reference_total)) {
    perror("Error when unmapping the reference file");
    goto closeexit;
  }

  close(reference_fd);
  return 0;

 usage:
  fprintf(stderr,
	  "usage: gen-hash <reference-genome>\n"
	  "<reference-genome> is the name of a packed binary reference sequence.\n");
 closeexit:
  if (reference_fd >= 0) close(reference_fd);
  return 1;
}
