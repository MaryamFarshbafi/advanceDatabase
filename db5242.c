/*
  CSE 5242 Project 2, Fall 2024

  See class project handout for more extensive documentation.
*/

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>
#include <sys/ioctl.h>
#include <sys/time.h>
#include <asm/unistd.h>
#include <immintrin.h>

/* uncomment out the following DEBUG line for debug info, for experiment comment the DEBUG line  */
// #define DEBUG

int64_t low_bin_search(int64_t *data, int64_t size, int64_t target);
void bulk_bin_search(int64_t *data, int64_t size, int64_t *searchkeys, int64_t numsearches, int64_t *results, int repeats);
int64_t low_bin_nb_arithmetic(int64_t *data, int64_t size, int64_t target);
int64_t low_bin_nb_mask(int64_t *data, int64_t size, int64_t target);
void low_bin_nb_4x(int64_t *data, int64_t size, int64_t *targets, int64_t *right);
void low_bin_nb_simd(int64_t *data, int64_t size, __m256i target, __m256i *result);

// test commit
/* compare two int64_t values - for use with qsort */
static int compare(const void *p1, const void *p2)
{
  int a, b;
  a = *(int64_t *)p1;
  b = *(int64_t *)p2;
  if (a < b)
    return -1;
  if (a == b)
    return 0;
  return 1;
}

/* initialize searches and data - data is sorted and searches is a random permutation of data */
int init(int64_t *data, int64_t *searches, int count)
{
  for (int64_t i = 0; i < count; i++)
  {
    searches[i] = random();
    data[i] = searches[i] + 1;
  }
  qsort(data, count, sizeof(int64_t), compare);
}

/* initialize outer probes of band join */
int band_init(int64_t *outer, int64_t size)
{
  for (int64_t i = 0; i < size; i++)
  {
    outer[i] = random();
  }
}

inline int64_t simple_binary_search(int64_t *data, int64_t size, int64_t target)
{
  int64_t left = 0;
  int64_t right = size;
  int64_t mid;

  while (left <= right)
  {
    mid = (left + right) / 2; /* ignore possibility of overflow of left+right */
    if (data[mid] == target)
      return mid;
    if (data[mid] < target)
      left = mid + 1;
    else
      right = mid - 1;
  }
  return -1; /* no match */
}

inline int64_t low_bin_search(int64_t *data, int64_t size, int64_t target)
{
  /* this binary search variant
     (a) does only one comparison in the inner loop
     (b) doesn't require an exact match; instead it returns the index of the first key >= the search key.
         That's good in a DB context where we might be doing a range search, and using binary search to
   identify the first key in the range.
     (c) If the search key is bigger than all keys, it returns size.
  */
  int64_t left = 0;
  int64_t right = size;
  int64_t mid;

  while (left < right)
  {
    mid = (left + right) / 2; /* ignore possibility of overflow of left+right */
    if (data[mid] >= target)
      right = mid;
    else
      left = mid + 1;
  }
  return right;
}

inline int64_t low_bin_nb_arithmetic(int64_t *data, int64_t size, int64_t target)
{
  /* this binary search variant
     (a) does no comparisons in the inner loop by using multiplication and addition to convert control dependencies
         to data dependencies
     (b) doesn't require an exact match; instead it returns the index of the first key >= the search key.
         That's good in a DB context where we might be doing a range search, and using binary search to
   identify the first key in the range.
     (c) If the search key is bigger than all keys, it returns size.
  */
  int64_t left = 0;
  int64_t right = size;
  int64_t mid;

  while (left < right)
  {
    /* YOUR CODE HERE */
    mid = (left + (right - left)) / 2;
    int64_t condition = (data[mid] >= target);
    right = condition * mid + (1 - condition) * right;
    left = (1 - condition) * (mid + 1) + condition * left;
  }
  return right;
}

inline int64_t low_bin_nb_mask(int64_t *data, int64_t size, int64_t target)
{
  /* this binary search variant
     (a) does no comparisons in the inner loop by using bit masking operations to convert control dependencies
         to data dependencies
     (b) doesn't require an exact match; instead it returns the index of the first key >= the search key.
         That's good in a DB context where we might be doing a range search, and using binary search to
   identify the first key in the range.
     (c) If the search key is bigger than all keys, it returns size.
  */
  int64_t left = 0;
  int64_t right = size;
  /* YOUR CODE HERE */
  int64_t mid;

  while (left < right)
  {
    mid = (left + (right - left)) / 2;
    int64_t mask = -(data[mid] >= target);
    right = (mask & mid) | (~mask & right);
    left = (~mask & (mid + 1)) | (mask & left);

    return right;
  }
}

int64_t low_bin_nb_mask_call(int64_t *data, int64_t size, int64_t target)
{
  /* this binary search variant
     (a) does no comparisons in the inner loop by using bit masking operations to convert control dependencies
         to data dependencies
     (b) doesn't require an exact match; instead it returns the index of the first key >= the search key.
         That's good in a DB context where we might be doing a range search, and using binary search to
   identify the first key in the range.
     (c) If the search key is bigger than all keys, it returns size.
  */
  int64_t left = 0;
  int64_t right = size;
  /* YOUR CODE HERE */
  int64_t mid;

  while (left < right)
  {
    mid = (left + (right - left)) / 2;
    int64_t mask = -(data[mid] >= target);
    right = (mask & mid) | (~mask & right);
    left = (~mask & (mid + 1)) | (mask & left);

    return right;
  }
}

inline void low_bin_nb_4x(int64_t *data, int64_t size, int64_t *targets, int64_t *right)
{
/* this binary search variant                                                                                       \
   (a) does no comparisons in the inner loop by using bit masking operations instead                                \
   (b) doesn't require an exact match; instead it returns the index of the first key >= the search key.             \
       That's good in a DB context where we might be doing a range search, and using binary search to               \
 identify the first key in the range.                                                                               \
   (c) If the search key is bigger than all keys, it returns size.                                                  \
   (d) does 4 searches at the same time in an interleaved fashion, so that an out-of-order processor                \
       can make progress even when some instructions are still waiting for their operands to be ready.              \
                                                                                                                  \ \
   Note that we're using the result array as the "right" elements in the search so no need for a return statement   \
*/                                                                                                                  \
 int64_t left[4] = {0, 0, 0, 0};
  int64_t mid[4];
  memcpy(right, (int64_t[]){size, size, size, size}, 4 * sizeof(int64_t));

  while (left[0] < right[0] || left[1] < right[1] || left[2] < right[2] || left[3] < right[3])
  {
    for (int i = 0; i < 4; i++)
    {
      if (left[i] < right[i])
      {
        mid[i] = (left[i] + (right[i] - left[i]) / 2);
        int64_t mask = -(data[mid[i]] >= targets[i]);
        right[i] = (mask & mid[i]) | (~mask & right[i]);
        left[i] = (~mask & (mid[i] + 1)) | (mask & left[i]);
      }
    }
  }
}

/* The following union type is handy to output the contents of AVX512 data types */
union int8x4
{
  __m256i a;
  int64_t b[4];
};

void printavx(char *name, __m256i v)
{
  union int8x4 n;

  n.a = v;
  printf("Value in %s is [%ld %ld %ld %ld ]\n", name, n.b[0], n.b[1], n.b[2], n.b[3]);
}

/*
 * Optinal for using AVX-512

  union int8x8 {
    __m512i a;
    int64_t b[8];
  };

  void printavx512(char* name, __m512i v) {
    union int8x4 n;

    n.a=v;
    printf("Value in %s is [%ld %ld %ld %ld %ld %ld %ld %ld ]\n",name,n.b[0],n.b[1],n.b[2],n.b[3]);
  }

 */

 inline void low_bin_nb_simd(int64_t *data, int64_t size, __m256i target, __m256i *result)
{
  /* this binary search variant
     (a) does no comparisons in the inner loop by using bit masking operations instead
     (b) doesn't require an exact match; instead it returns the index of the first key >= the search key.
         That's good in a DB context where we might be doing a range search, and using binary search to
   identify the first key in the range.
     (c) If the search key is bigger than all keys, it returns size.
     (d) does 4 searches at the same time using AVX2 intrinsics

     See https://www.intel.com/content/www/us/en/docs/cpp-compiler/developer-guide-reference/2021-8/intrinsics-for-avx2.html
     for documentation of the AVX512 intrinsics

     Note that we're using the result array as the "right" elements in the search, and that searchkey is being passed
     as an __m256i value rather than via a pointer.
  */

  /* YOUR CODE HERE */
  __m256i left = _mm256_set1_epi64x(0);
  __m256i right = _mm256_set1_epi64x(size);
  __m256i one = _mm256_set1_epi64x(1);
  __m256i mid;
  while (_mm256_movemask_epi8(_mm256_cmpgt_epi64(right, left)) != 0)
  {
    mid = _mm256_add_epi64(left, _mm256_srli_epi64(_mm256_sub_epi64(right, left), 1));
    __m256i mid_values = _mm256_i64gather_epi64((long long int *)data, mid, 8);
    __m256i cmp = _mm256_cmpgt_epi64(target, mid_values);
    left = _mm256_blendv_epi8(left, _mm256_add_epi64(mid, one), cmp);
    right = _mm256_blendv_epi8(mid, right, cmp);

    // now it needs to check the condition terminatuon
    __m256i finished = _mm256_cmpeq_epi64(left, right);
    if (_mm256_testz_si256(finished, finished))
    {
      break;
    }
  }
  *result = right;
}

void bulk_bin_search(int64_t *data, int64_t size, int64_t *searchkeys, int64_t numsearches, int64_t *results, int repeats)
{
  for (int j = 0; j < repeats; j++)
  {
    /* Function to test a large number of binary searches

       we might need repeats>1 to make sure the events we're measuring are not dominated by various
       overheads, particularly for small values of size and/or numsearches

       we assume that we want exactly "size" searches, where "size" is the length if the searchkeys array
     */
    for (int64_t i = 0; i < numsearches; i++)
    {
#ifdef DEBUG
      printf("Searching for %ld...\n", searchkeys[i]);
#endif

      // Uncomment one of the following to measure it
      results[i] = low_bin_search(data, size, searchkeys[i]);
      // results[i] = low_bin_nb_arithmetic(data,size,searchkeys[i]);
      // results[i] = low_bin_nb_mask(data,size,searchkeys[i]);

#ifdef DEBUG
      printf("Result is %ld\n", results[i]);
#endif
    }
  }
}

void bulk_bin_search_4x(int64_t *data, int64_t size, int64_t *searchkeys, int64_t numsearches, int64_t *results, int repeats)
{
  register __m256i searchkey_4x;

  for (int j = 0; j < repeats; j++)
  {
    /* Function to test a large number of binary searches using one of the 8x routines

       we might need repeats>1 to make sure the events we're measuring are not dominated by various
       overheads, particularly for small values of size and/or numsearches

       we assume that we want exactly "size" searches, where "size" is the length if the searchkeys array
     */
    int64_t extras = numsearches % 4;
    for (int64_t i = 0; i < numsearches - extras; i += 4)
    {
#ifdef DEBUG
      printf("Searching for %ld %ld %ld %ld  ...\n",
             searchkeys[i], searchkeys[i + 1], searchkeys[i + 2], searchkeys[i + 3]);
#endif

      // Uncomment one of the following depending on which routine you want to profile

      // Algorithm A
      low_bin_nb_4x(data, size, &searchkeys[i], &results[i]);

      // Algorithm B
      // searchkey_4x = _mm256_loadu_si256((__m256i *)&searchkeys[i]);
      // low_bin_nb_simd(data,size,searchkey_4x,(__m256i *)&results[i]);

#ifdef DEBUG
      printf("Result is %ld %ld %ld %ld  ...\n",
             results[i], results[i + 1], results[i + 2], results[i + 3]);
#endif
    }
    /* a little bit more work if numsearches is not a multiple of 8 */
    for (int64_t i = numsearches - extras; i < numsearches; i++)
    {

      results[i] = low_bin_nb_mask(data, size, searchkeys[i]);
    }
  }
}

int64_t band_join(int64_t *inner, int64_t inner_size, int64_t *outer, int64_t outer_size, int64_t *inner_results, int64_t *outer_results, int64_t result_size, int64_t bound)
{
  /* In a band join we want matches within a range of values.  If p is the probe value from the outer table, then all
     reccords in the inner table with a key in the range [p-bound,p+bound] inclusive should be part of the result.

     Results are returned via two arrays. outer_results stores the index of the outer table row that matches, and
     inner_results stores the index of the inner table row that matches.  result_size tells you the size of the
     output array that has been allocated. You should make sure that you don't exceed this size.  If there are
     more results than can fit in the result arrays, then return early with just a prefix of the results in the result
     arrays. The return value of the function should be the number of output results.

  */

  /* YOUR CODE HERE */
  /*input array:inner represent the column salary
  input array: outter represent request column
  output array: inner result
  output array: output result
  -it should be Use low_bin_nb_4x to show 4 binary searches simultaneously to overlap delay.
  -Use low_bin_nb_mask for any remaining searches.
  */
  int64_t count = 0; // track the number of result will be added
  // process four outer element
  int64_t index = 0;
  while (index + 4 <= outer_size)
  {
    int64_t lower_bound[4], upper_bound[4];
    // preper target range for the bash
    int64_t target_lower[4];
    int64_t target_upper[4];
    for (int i = 0; i < 4; i++)
    {
      target_lower[i] = outer[index + i] - bound;
      target_upper[i] = outer[index + i] + bound;
    }
    // show binary search for lower and upper bound
    low_bin_nb_4x(inner, inner_size, target_lower, lower_bound);
    low_bin_nb_4x(inner, inner_size, target_upper, upper_bound);

    // collect results for the batch
    for (int i = 0; i < 4; i++){
      for (int64_t k = lower_bound[i]; k < upper_bound[i]; k++)
      {
        if (count >= result_size)
        {
          return count;
        }
        inner_results[count] = k;
        outer_results[count] = index + i;
        count++;
      }
    }
    index += 4;
  }
  // process ramining outer element
  for (int i = 0; i < outer_size; i++)
  {
    int64_t lower_bound = low_bin_nb_mask(inner, inner_size, outer[i] - bound);
    int64_t upper_bound = low_bin_nb_mask(inner, inner_size, outer[i] + bound);
    for (int64_t j = lower_bound; j < upper_bound; j++)
    {
      if (count >= result_size)
      {
        return count;
        inner_results[count] = j;
        inner_results[count] = i;
        count++;
      }
    }
  }

  return count;
}

int64_t band_join_simd(int64_t *inner, int64_t inner_size, int64_t *outer, int64_t outer_size, int64_t *inner_results, int64_t *outer_results, int64_t result_size, int64_t bound)
{
  /* In a band join we want matches within a range of values.  If p is the probe value from the outer table, then all
     reccords in the inner table with a key in the range [p-bound,p+bound] inclusive should be part of the result.

     Results are returned via two arrays. outer_results stores the index of the outer table row that matches, and
     inner_results stores the index of the inner table row that matches.  result_size tells you the size of the
     output array that has been allocated. You should make sure that you don't exceed this size.  If there are
     more results than can fit in the result arrays, then return early with just a prefix of the results in the result
     arrays. The return value of the function should be the number of output results.

     To do the binary search, you could use the low_bin_nb_simd you just implemented to search for the lower bounds in parallel

     Once you've found the lower bounds, do the following for each of the 4 search keys in turn:
        scan along the sorted inner array, generating outputs for each match, and making sure not to exceed the output array bounds.

     This inner scanning code does not have to use SIMD.
  */
  int64_t output_count = 0;

  for (int64_t outer_idx = 0; outer_idx < outer_size; outer_idx += 4) {
      // Load the next 4 elements from the outer array
      int64_t keys[4] = {0};
      for (int i = 0; i < 4 && outer_idx + i < outer_size; ++i) {
          keys[i] = outer[outer_idx + i];
      }

      // Define the lower and upper bounds
      __m256i target_low = _mm256_sub_epi64(_mm256_loadu_si256((__m256i *)keys), _mm256_set1_epi64x(bound));
      __m256i target_high = _mm256_add_epi64(_mm256_loadu_si256((__m256i *)keys), _mm256_set1_epi64x(bound));

      // Perform binary search for the lower and upper bounds
      __m256i lower_indices, upper_indices;
      low_bin_nb_simd(inner, inner_size, target_low, &lower_indices);
      low_bin_nb_simd(inner, inner_size, target_high, &upper_indices);

      // Iterate over each key
      for (int i = 0; i < 4 && outer_idx + i < outer_size; ++i) {
          int64_t lower = ((int64_t *)&lower_indices)[i];
          int64_t upper = ((int64_t *)&upper_indices)[i];

          // Scan the inner array within the bounds
          for (int64_t inner_idx = lower; inner_idx < inner_size && inner_idx <= upper; ++inner_idx) {
              if (inner[inner_idx] >= keys[i] - bound && inner[inner_idx] <= keys[i] + bound) {
                  // Add results to the output arrays
                  if (output_count < result_size) {
                      inner_results[output_count] = inner_idx;
                      outer_results[output_count] = outer_idx + i;
                      output_count++;
                  } else {
                      return output_count; // Return early if the result buffer is full
                  }
              }
          }
      }
  }
  return output_count;
}

void write_csv_header(FILE *file) {
    fprintf(file, "Function,Total Time (microseconds),Avg Time Per Search (microseconds),Total Results,Avg Matches Per Output Record\n");
}

void write_to_csv(FILE *file, const char *function_name, long total_time, double avg_time_per_search, long total_results, double avg_matches_per_record) {
    fprintf(file, "%s,%ld,%f,%ld,%f\n", function_name, total_time, avg_time_per_search, total_results, avg_matches_per_record);
}

int main(int argc, char *argv[])
{
  long long counter;
  int64_t arraysize, outer_size, result_size;
  int64_t bound;
  int64_t *data, *queries, *results;
  int ret;
  struct timeval before, after;
  int repeats;
  int64_t total_results;

  // band-join arrays
  int64_t *outer, *outer_results, *inner_results;

  if (argc >= 5)
  {
    arraysize = atoi(argv[1]);
    outer_size = atoi(argv[2]);
    result_size = atoi(argv[3]);
    bound = atoi(argv[4]);
  }
  else
  {
    fprintf(stderr, "Usage: db5242 inner_size outer_size result_size bound <repeats>\n");
    exit(EXIT_FAILURE);
  }

  if (argc >= 6)
    repeats = atoi(argv[5]);
  else
  {
    repeats = 1;
  }

  /* allocate the array and the queries for searching */
  ret = posix_memalign((void **)&data, 64, arraysize * sizeof(int64_t));
  if (ret)
  {
    fprintf(stderr, "Memory allocation error.\n");
    exit(EXIT_FAILURE);
  }
  ret = posix_memalign((void **)&queries, 64, arraysize * sizeof(int64_t));
  if (ret)
  {
    fprintf(stderr, "Memory allocation error.\n");
    exit(EXIT_FAILURE);
  }
  ret = posix_memalign((void **)&results, 64, arraysize * sizeof(int64_t));
  if (ret)
  {
    fprintf(stderr, "Memory allocation error.\n");
    exit(EXIT_FAILURE);
  }

  /* allocate the outer array and output arrays for band-join */
  ret = posix_memalign((void **)&outer, 64, outer_size * sizeof(int64_t));
  if (ret)
  {
    fprintf(stderr, "Memory allocation error.\n");
    exit(EXIT_FAILURE);
  }
  ret = posix_memalign((void **)&outer_results, 64, result_size * sizeof(int64_t));
  if (ret)
  {
    fprintf(stderr, "Memory allocation error.\n");
    exit(EXIT_FAILURE);
  }
  ret = posix_memalign((void **)&inner_results, 64, result_size * sizeof(int64_t));
  if (ret)
  {
    fprintf(stderr, "Memory allocation error.\n");
    exit(EXIT_FAILURE);
  }

  /* code to initialize data structures goes here so that it is not included in the timing measurement */
  init(data, queries, arraysize);
  band_init(outer, outer_size);

#ifdef DEBUG
  /* show the arrays */
  printf("data: ");
  for (int64_t i = 0; i < arraysize; i++)
    printf("%ld ", data[i]);
  printf("\n");
  printf("queries: ");
  for (int64_t i = 0; i < arraysize; i++)
    printf("%ld ", queries[i]);
  printf("\n");
  printf("outer: ");
  for (int64_t i = 0; i < outer_size; i++)
    printf("%ld ", outer[i]);
  printf("\n");
#endif

  /* now measure... */
  FILE *file = fopen("performance_results.csv", "w");
  if (!file) {
      printf("Failed to open file for writing\n");
      return 1;
  }
  
  write_csv_header(file);

  // Measure and log for bulk_bin_search
  gettimeofday(&before, NULL);
  bulk_bin_search(data, arraysize, queries, arraysize, results, repeats);
  gettimeofday(&after, NULL);

  long bulk_bin_search_time = (after.tv_sec - before.tv_sec) * 1000000 + (after.tv_usec - before.tv_usec);
  double bulk_bin_search_avg = 1.0 * bulk_bin_search_time / arraysize / repeats;
  write_to_csv(file, "bulk_bin_search", bulk_bin_search_time, bulk_bin_search_avg, 0, 0);

  // Measure and log for bulk_bin_search_4x
  gettimeofday(&before, NULL);
  bulk_bin_search_4x(data, arraysize, queries, arraysize, results, repeats);
  gettimeofday(&after, NULL);

  long bulk_bin_search_4x_time = (after.tv_sec - before.tv_sec) * 1000000 + (after.tv_usec - before.tv_usec);
  double bulk_bin_search_4x_avg = 1.0 * bulk_bin_search_4x_time / arraysize / repeats;
  write_to_csv(file, "bulk_bin_search_4x", bulk_bin_search_4x_time, bulk_bin_search_4x_avg, 0, 0);

  // Measure and log for band_join
  gettimeofday(&before, NULL);
  total_results = band_join(data, arraysize, outer, outer_size, inner_results, outer_results, result_size, bound);
  gettimeofday(&after, NULL);

  long band_join_time = (after.tv_sec - before.tv_sec) * 1000000 + (after.tv_usec - before.tv_usec);
  double band_join_avg_matches = 1.0 * total_results / (1.0 + outer_results[total_results - 1]);
  write_to_csv(file, "band_join", band_join_time, 0, total_results, band_join_avg_matches);

  fclose(file);
  printf("Finished writing to performance_results.csv\n");
  return 0;

#ifdef DEBUG
  /* show the band_join results */
  printf("band_join results: ");
  for (int64_t i = 0; i < total_results; i++)
    printf("(%ld,%ld) ", outer_results[i], inner_results[i]);
  printf("\n");
#endif
}
