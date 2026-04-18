/*
 * Standalone reproducer of NCBI ALIGN_EX (blast_gapalign.c:374).
 * Strips the NCBI struct dependencies (GapPrelimEditBlock, BlastGapAlignStruct,
 * BlastScoreBlk, GapStateArrayStruct) down to what the DP recurrence actually
 * touches, so we can instrument it with tracehash-rs's C helpers and diff
 * row-by-row against Rust's `align_ex` in `src/traceback.rs`.
 *
 * The DP logic is a line-by-line copy of NCBI's ALIGN_EX. No tie-break rules
 * or predecessor order have been modified.
 *
 * Build:
 *   cc -O0 -g -Wall -Wextra \
 *      -I /home/mahogny/.cargo/registry/src/index.crates.io-XXX/tracehash-rs-0.1.0/c \
 *      ncbi_align_ex.c \
 *      /home/mahogny/.cargo/registry/src/index.crates.io-XXX/tracehash-rs-0.1.0/c/tracehash_c.c \
 *      -o ncbi_align_ex
 *
 * Run:
 *   TRACEHASH_OUT=/tmp/ncbi-trace.tsv TRACEHASH_SIDE=ncbi \
 *   TRACEHASH_RUN_ID=adj_del_ins ./ncbi_align_ex
 */
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tracehash_c.h"

/* NCBI typedefs (matching core/blast_def.h). */
typedef unsigned char Uint1;
typedef int Int4;

#define INT4_MIN_HALF ((Int4)(-1073741824)) /* i32::MIN / 2, matches Rust MININT */
#define FENCE_SENTRY 0x7f /* unused here since we emit no fence */

/* Script/op constants (blast_gapalign.c:363-371, gapinfo.h:45-51). */
#define SCRIPT_SUB 3
#define SCRIPT_GAP_IN_A 0
#define SCRIPT_GAP_IN_B 6
#define SCRIPT_OP_MASK 0x07
#define SCRIPT_EXTEND_GAP_A 0x10
#define SCRIPT_EXTEND_GAP_B 0x40

typedef struct {
  Int4 best;
  Int4 best_gap;
} BlastGapDP;

#define MAX_INT(a, b) ((a) > (b) ? (a) : (b))

/* Stripped ALIGN_EX. Same DP, no traceback assembly (we only need DP state). */
static Int4 align_ex(const Uint1 *A, const Uint1 *B, Int4 M, Int4 N,
                     Int4 *a_offset, Int4 *b_offset, const Int4 *matrix,
                     size_t matrix_stride, Int4 gap_open, Int4 gap_extend,
                     Int4 x_dropoff, int reverse_sequence) {
  Int4 i;
  Int4 a_index;
  Int4 b_index, b_size, first_b_index, last_b_index, b_increment;
  const Uint1 *b_ptr;
  BlastGapDP *score_array;
  Int4 gap_open_extend;
  Int4 best_score;
  const Int4 *matrix_row = NULL;
  Int4 score, score_gap_row, score_gap_col, next_score;
  Uint1 script, script_col, script_row;
  Uint1 **edit_script;
  Int4 *edit_start_offset;
  Uint1 *edit_script_row;
  Int4 edit_script_num_rows = 100;
  Int4 num_extra_cells;

  gap_open_extend = gap_open + gap_extend;
  if (x_dropoff < gap_open_extend)
    x_dropoff = gap_open_extend;
  if (N <= 0 || M <= 0)
    return 0;

  if (gap_extend > 0)
    num_extra_cells = x_dropoff / gap_extend + 3;
  else
    num_extra_cells = N + 3;

  Int4 alloc = N + num_extra_cells + 10;
  score_array = (BlastGapDP *)calloc(alloc, sizeof(BlastGapDP));
  edit_script = (Uint1 **)malloc(sizeof(Uint1 *) * edit_script_num_rows);
  edit_start_offset = (Int4 *)malloc(sizeof(Int4) * edit_script_num_rows);
  edit_script[0] = (Uint1 *)calloc(N + num_extra_cells + 10, 1);
  edit_start_offset[0] = 0;
  edit_script_row = edit_script[0];

  score = -gap_open_extend;
  score_array[0].best = 0;
  score_array[0].best_gap = -gap_open_extend;
  for (i = 1; i <= N; i++) {
    if (score < -x_dropoff)
      break;
    score_array[i].best = score;
    score_array[i].best_gap = score - gap_open_extend;
    score -= gap_extend;
    edit_script_row[i] = SCRIPT_GAP_IN_A;
  }
  b_size = i;
  best_score = 0;
  first_b_index = 0;
  b_increment = reverse_sequence ? -1 : 1;
  *a_offset = 0;
  *b_offset = 0;

  for (a_index = 1; a_index <= M; a_index++) {
    if (a_index >= edit_script_num_rows) {
      edit_script_num_rows *= 2;
      edit_script = (Uint1 **)realloc(edit_script,
                                      edit_script_num_rows * sizeof(Uint1 *));
      edit_start_offset =
          (Int4 *)realloc(edit_start_offset, edit_script_num_rows * sizeof(Int4));
    }
    edit_script[a_index] =
        (Uint1 *)calloc(N + num_extra_cells + 10 - first_b_index + 1, 1);
    edit_start_offset[a_index] = first_b_index;
    edit_script_row = edit_script[a_index] - first_b_index;

    if (reverse_sequence)
      matrix_row = &matrix[A[M - a_index] * matrix_stride];
    else
      matrix_row = &matrix[A[a_index] * matrix_stride];

    if (reverse_sequence)
      b_ptr = &B[N - first_b_index];
    else
      b_ptr = &B[first_b_index];

    score = INT4_MIN_HALF;
    score_gap_row = INT4_MIN_HALF;
    last_b_index = first_b_index;

    for (b_index = first_b_index; b_index < b_size; b_index++) {
      b_ptr += b_increment;
      score_gap_col = score_array[b_index].best_gap;
      Int4 b_letter = *b_ptr;
      next_score = score_array[b_index].best + matrix_row[b_letter];

      script = SCRIPT_SUB;
      script_col = SCRIPT_EXTEND_GAP_B;
      script_row = SCRIPT_EXTEND_GAP_A;

      if (score < score_gap_col) {
        script = SCRIPT_GAP_IN_B;
        score = score_gap_col;
      }
      if (score < score_gap_row) {
        script = SCRIPT_GAP_IN_A;
        score = score_gap_row;
      }

      if (best_score - score > x_dropoff) {
        if (first_b_index == b_index)
          first_b_index++;
        else
          score_array[b_index].best = INT4_MIN_HALF;
      } else {
        last_b_index = b_index;
        if (score > best_score) {
          best_score = score;
          *a_offset = a_index;
          *b_offset = b_index;
        }

        score_gap_row -= gap_extend;
        score_gap_col -= gap_extend;
        if (score_gap_col < (score - gap_open_extend)) {
          score_array[b_index].best_gap = score - gap_open_extend;
        } else {
          score_array[b_index].best_gap = score_gap_col;
          script += script_col;
        }

        if (score_gap_row < (score - gap_open_extend))
          score_gap_row = score - gap_open_extend;
        else
          script += script_row;

        score_array[b_index].best = score;
      }

      score = next_score;
      edit_script_row[b_index] = script;
    }

    /* Emit one tracehash event per row, matching Rust's `align_ex_row`. */
    {
      TraceHashCall call = tracehash_begin("align_ex_row", __FILE__, __LINE__);
      tracehash_input_u64(&call, (uint64_t)a_index);
      tracehash_input_u64(&call, (uint64_t)first_b_index);
      tracehash_input_u64(&call, (uint64_t)last_b_index);
      tracehash_input_u64(&call, (uint64_t)b_size);
      tracehash_input_i64(&call, (int64_t)best_score);
      tracehash_input_i64(&call, (int64_t)*a_offset);
      tracehash_input_i64(&call, (int64_t)*b_offset);
      /* Pack score_array[0..b_size] as interleaved i32 little-endian bytes.
         Pre-row value = DP state we just computed for this row. */
      size_t buf_len = (size_t)b_size * 8;
      uint8_t *buf = (uint8_t *)malloc(buf_len);
      for (Int4 bi = 0; bi < b_size; bi++) {
        int32_t vb = (int32_t)score_array[bi].best;
        int32_t vg = (int32_t)score_array[bi].best_gap;
        memcpy(buf + (size_t)bi * 8 + 0, &vb, 4);
        memcpy(buf + (size_t)bi * 8 + 4, &vg, 4);
      }
      tracehash_input_bytes(&call, buf, buf_len);
      free(buf);
      tracehash_input_bytes(&call, edit_script_row + first_b_index,
                            (size_t)(b_size - first_b_index));
      tracehash_finish(&call);
    }

    if (first_b_index == b_size)
      break;

    if (last_b_index < b_size - 1) {
      b_size = last_b_index + 1;
    } else {
      while (score_gap_row >= (best_score - x_dropoff) && b_size <= N) {
        score_array[b_size].best = score_gap_row;
        score_array[b_size].best_gap = score_gap_row - gap_open_extend;
        score_gap_row -= gap_extend;
        edit_script_row[b_size] = SCRIPT_GAP_IN_A;
        b_size++;
      }
    }

    if (b_size <= N) {
      score_array[b_size].best = INT4_MIN_HALF;
      score_array[b_size].best_gap = INT4_MIN_HALF;
      b_size++;
    }
  }

  /* Traceback walk — copy of blast_gapalign.c:682-727. Emit one tracehash
   * event per op so the Rust side can match exactly. */
  a_index = *a_offset;
  b_index = *b_offset;
  script = SCRIPT_SUB;
  Int4 step = 0;
  while (a_index > 0 || b_index > 0) {
    Uint1 next_script =
        edit_script[a_index][b_index - edit_start_offset[a_index]];

    switch (script) {
    case SCRIPT_GAP_IN_A:
      script = next_script & SCRIPT_OP_MASK;
      if (next_script & SCRIPT_EXTEND_GAP_A)
        script = SCRIPT_GAP_IN_A;
      break;
    case SCRIPT_GAP_IN_B:
      script = next_script & SCRIPT_OP_MASK;
      if (next_script & SCRIPT_EXTEND_GAP_B)
        script = SCRIPT_GAP_IN_B;
      break;
    default:
      script = next_script & SCRIPT_OP_MASK;
      break;
    }

    {
      TraceHashCall tcall = tracehash_begin("align_ex_tb", __FILE__, __LINE__);
      tracehash_input_u64(&tcall, (uint64_t)step);
      tracehash_input_u64(&tcall, (uint64_t)a_index);
      tracehash_input_u64(&tcall, (uint64_t)b_index);
      tracehash_input_u64(&tcall, (uint64_t)script);
      tracehash_finish(&tcall);
    }
    step++;

    if (script == SCRIPT_GAP_IN_A) {
      b_index--;
    } else if (script == SCRIPT_GAP_IN_B) {
      a_index--;
    } else {
      a_index--;
      b_index--;
    }
  }

  for (Int4 r = 0; r < (Int4)edit_script_num_rows; r++) {
    if (r <= *a_offset && edit_script[r])
      free(edit_script[r]);
  }
  free(edit_script);
  free(edit_start_offset);
  free(score_array);
  return best_score;
}

/* Build a BLASTNA match/mismatch matrix for reward=1, penalty=-3.
 * BLASTNA indices 0-3 = A,C,G,T. 4-15 = ambiguity codes. We only need the
 * first 16 rows × 16 cols for this parity probe (no ambiguity in the test
 * sequences). */
static void build_blastna_matrix(Int4 *matrix, size_t stride, Int4 reward,
                                 Int4 penalty) {
  for (size_t i = 0; i < 16; i++)
    for (size_t j = 0; j < 16; j++)
      matrix[i * stride + j] = (i < 4 && j < 4 && i == j) ? reward : penalty;
}

static Uint1 encode(char c) {
  switch (c) {
  case 'A':
    return 0;
  case 'C':
    return 1;
  case 'G':
    return 2;
  case 'T':
    return 3;
  default:
    return 15;
  }
}

static void encode_all(const char *s, Uint1 *out) {
  for (size_t i = 0; s[i]; i++)
    out[i] = encode(s[i]);
}

int main(void) {
  /* Sequences from tests/integration.rs adjacent_del_ins case. */
  const char *q_str =
      "ACGTTGCAACGATCGTACGATTCGAGCTTAGGCTAGGGTAATCGGATCCTAGCTAGGCTAATCGATCGTAGCTAGCATCGAT";
  const char *s_str =
      "ACGTTGCAACGATCGTACGATTCGAGCTTAGGCTAATCGGATCCTAGCTAGGCTAATCGATCGTAGCTAGCATCGAT";
  size_t qlen = strlen(q_str);
  size_t slen = strlen(s_str);
  Uint1 *q = (Uint1 *)malloc(qlen + 1);
  Uint1 *s = (Uint1 *)malloc(slen + 1);
  encode_all(q_str, q);
  encode_all(s_str, s);

  /* Seed and scoring match the Rust reproducer exactly. */
  Int4 seed_q = 17;
  Int4 seed_s = 17;
  Int4 reward = 1, penalty = -3;
  Int4 gap_open = 5, gap_extend = 2;
  Int4 x_dropoff = 16;

  Int4 matrix[16 * 16];
  build_blastna_matrix(matrix, 16, reward, penalty);

  /* Left extension: reverse_sequence=1, M=seed_q+1, slices = query[..seed_q+1],
   * subject[..seed_s+1]. */
  Int4 left_a_off = 0, left_b_off = 0;
  Int4 score_left = align_ex(q, s, seed_q + 1, seed_s + 1, &left_a_off,
                             &left_b_off, matrix, 16, gap_open, gap_extend,
                             x_dropoff, /*reverse_sequence*/ 1);

  /* Right extension: reverse_sequence=0, M=qlen-seed_q-1, slices =
   * query[seed_q..], subject[seed_s..]. */
  Int4 right_a_off = 0, right_b_off = 0;
  Int4 score_right =
      align_ex(q + seed_q, s + seed_s, (Int4)qlen - seed_q - 1,
               (Int4)slen - seed_s - 1, &right_a_off, &right_b_off, matrix,
               16, gap_open, gap_extend, x_dropoff, /*reverse_sequence*/ 0);

  fprintf(stderr,
          "left: score=%d a_off=%d b_off=%d; right: score=%d a_off=%d b_off=%d; "
          "total=%d\n",
          score_left, left_a_off, left_b_off, score_right, right_a_off,
          right_b_off, score_left + score_right);

  free(q);
  free(s);
  return 0;
}
