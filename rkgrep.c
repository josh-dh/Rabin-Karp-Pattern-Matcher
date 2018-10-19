/* NYU Computer Systems Organization Lab 2
 * Rabin-Karp Substring Matching
*/
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <getopt.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <ctype.h>
#include <math.h>

#include "rkgrep.h"
#include "bloom.h"

#define PRIME 961748941

// calculate modulo addition, i.e. (a+b) % PRIME
long long
madd(long long a, long long b)
{
	return (a + b) % PRIME;
}

// calculate modulo substraction, i.e. (a-b) % PRIME
long long
msub(long long a, long long b)
{
	return (a>b)?(a-b):(a+PRIME-b);
}

// calculate modulo multiplication, i.e. (a*b) % PRIME
long long
mmul(long long a, long long b)
{
	return (a*b) % PRIME;
}

/* naive_substring_match returns number of positions in the document where
 * the pattern has been found.  In addition, it stores the first position 
 * where the pattern is found in the variable pointed to by first_match_ind
 *
 * Both its arguments "pattern" and "doc" are null-terminated C strings.
 */
int
naive_substring_match(const char *pattern, const char *doc, int *first_match_ind)
{
	int match_count = 0;
	int m = strlen(pattern);
	int doc_len = strlen(doc);

	for (int i = 0; i < doc_len - m; i++) { //test all len(doc) - len(pattern) indices
		int equals = 1;
		for (int j = 0; j < m; j++) {
			if (pattern[j] != doc[i + j]) { //not matching
				equals = 0;
				break;
			}
		}
		if (equals) {
			if (match_count == 0) {
				*first_match_ind = i; //if first match, set first_match_ind to i
			}
			match_count += 1;
		}
	}
	return match_count;
}

/* initialize the Rabin-karp hash computation by calculating 
 * and returning the RK hash over a charbuf of m characters, 
 * i.e. The return value should be 
 * 256^(m-1)*charbuf[0] + 256^(m-2)*charbuf[1] + ... + charbuf[m-1],
 * where 256^(m-1) means 256 raised to the power m-1.
*/
long long
rkhash_init(const char *charbuf, int m, long long *h)
{
	*h = 1;
	long long sum = 0;
	for (int i = m - 1; i >= 0; i--) {
		sum = madd(mmul(*h, charbuf[i]), sum);
		*h = mmul(*h, 256);
	}
	return sum; //sum of 256 powers times their corresponding characters
}


/* Given the rabin-karp hash value (curr_hash) over substring Y[i],Y[i+1],...,Y[i+m-1]
 * calculate the hash value over Y[i+1],Y[i+2],...,Y[i+m] = curr_hash * 256 - leftmost * h + rightmost
 * where h is 256 raised to the power m (and given as an argument).  
 */
long long 
rkhash_next(long long curr_hash, long long h, char leftmost, char rightmost)
{
	return madd(msub(mmul(curr_hash, 256), mmul(leftmost, h)), rightmost);
}

/* rk_substring_match returns the number of positions in the document "doc" where
 * the "pattern" has been found, using the Rabin-karp substring matching algorithm.
 * Both pattern and doc are null-terminated C strings. The function also stores
 * the first position where pattern is found in the int variable pointed to by first_match_ind
*/
int
rk_substring_match(const char *pattern, const char *doc, int *first_match_ind)
{
	//init numbers
	int m = strlen(pattern);
	long long h;
	long long h_decoy; //used when calling second rkhash_init to avoid modifying h again
	int sum_matches = 0;

	//init hash patterns
    long long curr_hash = rkhash_init(doc, m, &h);
    long long pattern_hash = rkhash_init(pattern, m, &h_decoy);

    //check hashes and verify matches
    for (int i = 0; i < strlen(doc) - m; i++) { //loop from i = 0 to i = strlen(doc) - m

    	if (curr_hash == pattern_hash) { //verify
    		int ver = 1;
    		for (int j = i; j < i + m; j++) {
    			if (pattern[j - i] != doc[j]) { //if one character does not match, set ver to false and break
    				ver = 0;
    				break;
    			}
    		}
    		if (ver) {
    			sum_matches++;
    		}
    		if (ver == 1) {
    			*first_match_ind = i;
    		}
    	}

    	//set curr_hash to next value
    	curr_hash = rkhash_next(curr_hash, h, doc[i], doc[i + m]);
    }

	return sum_matches;
}


/* rk_create_doc_bloom returns a pointer to a newly created bloom_filter. 
 * The new bloom filter is populated with all n-m+1 rabin-karp hashes for 
 * all the substrings of length m in "doc".
 */
bloom_filter *
rk_create_doc_bloom(int m, const char *doc, int bloom_size)
{
	//init filter
	bloom_filter *bf = bloom_init(bloom_size);
	//create hashes and add to filter
	long long h;
	long long curr_hash = rkhash_init(doc, m, &h);
	for (int i = 0; i < strlen(doc) - m + 1; i++) { //iterate hashes
		//add to filter
		bloom_add(bf, curr_hash);
		//advance hash
		curr_hash = rkhash_next(curr_hash, h, doc[i], doc[i + m]);
	}

	return bf;
}

/* rk_substring_match_using_bloom returns the total number of positions where "pattern" 
 * is found in "doc".  It performs the matching by first checking against the 
 * pre-populated bloom filter "bf" (which has been created by rk_create_doc_bloom on "doc")
 * If the pattern is not found in "bf", then the function immediately returns 0.
 * Otherwise, the function invokes rk_substring_match() function to find "pattern" in "doc".
*/
int
rk_substring_match_using_bloom(const char *pattern, const char *doc, bloom_filter *bf, int *first_match_ind)
{
	//init pattern hash
	long long h;
	long long pattern_hash = rkhash_init(pattern, strlen(pattern), &h);
	//query and return 0 if false
    if (bloom_query(bf, pattern_hash)) {
    	//substring match
    	rk_substring_match(pattern, doc, first_match_ind);
    } else {
    	return 0;
    }
}
