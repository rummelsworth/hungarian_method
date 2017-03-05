/*
 * Copyright 2010 William Rummler (w.a.rummler@gmail.com)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * This file contains a C implementation of the Hungarian method as described
 * in Chapter 11 of "Combinatorial Optimization: Algorithms and Complexity" by
 * Papadimitriou and Steiglitz. All pages and figures refer to the 1998 Dover
 * edition of the book, corrected and up-to-date with the 8 October 2000 errata
 * file located at <http://www.cs.princeton.edu/~ken/latest.pdf> as of 21
 * November 2010.
 */
#include "hungarian_method.h"
#include <stdlib.h>
#include <limits.h>

/*
 * Solely for the tracing function hm_print() defined below.
 */
#include <stdio.h>

/*
 * This implementation uses zero-based indexing with V={0,...,n-1} and
 * U={n,...,2n-1}, as opposed to the book where use of zero as a "blank"
 * value precludes zero-based indexing. Thus, we use a negative int as a
 * "blank" value.
 */
enum { blank = -1 };

/*
 * Convenient macros to reduce clutter, improve readability, and facilitate
 * translation of vertex labels to zero-based indices where needed.
 */
#define Q              ( hm->q )
#define A              ( hm->a )
#define N              ( hm->n )
#define V( i_ )        ( ( i_ ) )
#define U( j_ )        ( ( j_ ) - N )
#define MATE( i_ )     ( hm->mate[ V( i_ ) ] )
#define C( i_, j_ )    ( hm->c[ V( i_ ) * N + U( j_ ) ] )
#define ALPHA( i_ )    ( hm->alpha[ V( i_ ) ] )
#define BETA( j_ )     ( hm->beta[ U( j_ ) ] )
#define SLACK( j_ )    ( hm->slack[ U( j_ ) ] )
#define NHBOR( j_ )    ( hm->nhbor[ U( j_ ) ] )
#define COUNT( i_ )    ( hm->count[ V( i_ ) ] )
#define EXPOSED( i_ )  ( hm->exposed[ V( i_ ) ] )
#define LABEL( i_ )    ( hm->label[ V( i_ ) ] )
#define EACH_V( i_ )   ( i_ = 0; i_ < N; ++i_ )
#define EACH_U( j_ )   ( j_ = N; j_ < 2 * N; ++j_ )

/*
 * Our basic data structures.
 */
typedef struct stack_    stack;
typedef struct arc_      arc;
typedef struct arc_list_ arc_list;
typedef struct hm_data_  hm_data;

#ifndef __cplusplus
/*
 * For convenience and clarity.
 */
typedef enum bool_ bool;
enum bool_ { false, true };
#endif

/*
 * The data structure corresponding to Q (see Figure 11-2).
 */
struct stack_
{
  int *data;
  int size;
};

static void
stack_push( stack *s, int datum )
{
  s->data[ s->size++ ] = datum;
}

static int
stack_pop( stack *s )
{
  return s->data[ --s->size ];
}

/*
 * The data structures corresponding to A (see Figure 11-2).
 */
struct arc_
{
  int x;
  int y;
};

struct arc_list_
{
  arc *data;
  int size;
};

static void
add_arc( arc_list *a, int x, int y )
{
  arc *pa = a->data + a->size++;
  pa->x = x;
  pa->y = y;
}

/*
 * This structure type holds all data pertinent to P&S's Hungarian method.
 * See Figure 11-2.
 */
struct hm_data_
{
  /*
   * Allocated and/or defined by the caller of hungarian_method().
   */
  int *mate;
  int *c;
  int n;
  /*
   * Allocated internally.
   */
  stack q;
  arc_list a;
  int *alpha;
  int *beta;
  int *slack;
  int *nhbor;
  int *count;
  int *exposed;
  int *label;
};

static void
hm_data_internal_free( hm_data *hm )
{
  if ( hm->q.data ) free( hm->q.data );
  if ( hm->a.data ) free( hm->a.data );
  if ( hm->alpha ) free( hm->alpha );
  if ( hm->beta ) free( hm->beta );
  if ( hm->slack ) free( hm->slack );
  if ( hm->nhbor ) free( hm->nhbor );
  if ( hm->count ) free( hm->count );
  if ( hm->exposed ) free( hm->exposed );
  if ( hm->label ) free( hm->label );
}

static hm_data *
hm_data_internal_allocate( hm_data *hm, int n )
{
  hm->q.data = NULL;
  hm->a.data = NULL;
  hm->alpha = NULL;
  hm->beta = NULL;
  hm->slack = NULL;
  hm->nhbor = NULL;
  hm->count = NULL;
  hm->exposed = NULL;
  hm->label = NULL;
  if ( ( hm->q.data = malloc( n * sizeof( *hm->q.data ) ) )
       && ( hm->a.data = malloc( n * n * sizeof( *hm->a.data ) ) )
       && ( hm->alpha = malloc( n * sizeof( *hm->alpha ) ) )
       && ( hm->beta = malloc( n * sizeof( *hm->beta ) ) )
       && ( hm->slack = malloc( n * sizeof( *hm->slack ) ) )
       && ( hm->nhbor = malloc( n * sizeof( *hm->nhbor ) ) )
       && ( hm->count = malloc( n * sizeof( *hm->count ) ) )
       && ( hm->exposed = malloc( n * sizeof( *hm->exposed ) ) )
       && ( hm->label = malloc( n * sizeof( *hm->label ) ) ) )
    {
      hm->n = n;
      return hm;
    }
  hm_data_internal_free( hm );
  return NULL;
}

/*
 * This function is for debugging purposes. It prints the algorithm's internal
 * state in a format similar to that of Example 11.1 (The matrix form of the
 * Hungarian method) beginning on page 252.
 *
 * The formatted output coded here is intended for small numbers.
 */
static void
hm_print( hm_data *hm )
{
  int i, j, k;
  printf( "\n a\\b |" );
  for EACH_U( j )
    {
      printf( "%3d ", BETA( j ) );
    }
  printf( "mate exposed label\n" );
  printf( "-----+" );
  for EACH_U( j )
    {
      printf( "----" );
    }
  printf( "------------------\n" );
  for EACH_V( i )
    {
      printf( "     |\n %3d |", ALPHA( i ) );
      for EACH_U( j )
        {
          printf( "%3d ", C( i, j ) );
        }
      printf( "%4d %7d %5d\n", MATE( i ), EXPOSED( i ), LABEL( i ) );
    }
  printf( "\nslack" );
  for EACH_U( j )
    {
      printf( " %3d", SLACK( j ) == INT_MAX ? -1 : SLACK( j ) );
    }
  printf( "\nnhbor" );
  for EACH_U( j )
    {
      printf( " %3d", NHBOR( j ) );
    }
  printf( "\n\nA = { " );
  for ( k = 0; k < A.size; ++k )
    {
      printf( "(%d,%d) ", A.data[ k ].x, A.data[ k ].y );
    }
  printf( "}\nQ = { " );
  for ( k = 0; k < Q.size; ++k )
    {
      printf( "%d ", Q.data[ k ] );
    }
  printf( "}\n\n" );
}

/*
 * See Figure 10-3, "The bipartite matching algorithm", page 224.
 *
 * Corresponds to "procedure augment(v)",
 * but is iterative instead of recursive.
 */
static void
hm_augment( hm_data *hm, int v )
{
  while ( LABEL( v ) != blank )
    {
      EXPOSED( LABEL( v ) ) = MATE( v );
      MATE( v ) = EXPOSED( v );
      MATE( EXPOSED( v ) ) = v;
      v = LABEL( v );
    }
  MATE( v ) = EXPOSED( v );
  MATE( EXPOSED( v ) ) = v;
}

/*
 * See Figure 11-2, "The Hungarian method", page 251.
 *
 * Corresponds to lines 7--8.
 */
static void
hm_initialize( hm_data *hm )
{
  int i, j;
  for EACH_V( i )
    {
      MATE( i ) = blank;
      ALPHA( i ) = 0;
    }
  for EACH_U( j )
    {
      MATE( j ) = blank;
      BETA( j ) = INT_MAX;
      for EACH_V( i )
        {
          if ( C( i, j ) < BETA( j ) )
            {
              BETA( j ) = C( i, j );
            }
        }
    }
}

/*
 * See Figure 11-2, "The Hungarian method", page 251.
 *
 * Corresponds to lines 12--17.
 */
static void
hm_construct_auxiliary_graph( hm_data *hm )
{
  int i, j;
  A.size = 0;
  for EACH_V( i )
    {
      EXPOSED( i ) = blank;
      LABEL( i ) = blank;
      /*
       * The following data structure is not included in the Figure 11-2
       * pseudo-code implementation. It has been added to account for
       * "labeling" on certain vertices described within Example 11.1 that
       * would otherwise be missing from the Figure 11-2 implementation.
       *
       * count[v] for any v \in V is equal to the size of the set
       * { u \in U : nhbor[u] = v }. When this set is non-empty, v is
       * considered to be "labeled". The use of this new data structure is
       * only to complete the conditional check on "labeled" statuses when
       * updating alpha within "procedure modify".
       */
      COUNT( i ) = 0;
    }
  for EACH_U( j )
    {
      SLACK( j ) = INT_MAX;
      /*
       * The following initialization of nhbor[] is necessary for proper usage
       * of the count[] array, whose addition and purpose is described above.
       */
      NHBOR( j ) = blank;
    }
  for EACH_V( i )
    {
      for EACH_U( j )
        {
          if ( ALPHA( i ) + BETA( j ) == C( i, j ) )
            {
              if ( MATE( j ) == blank )
                {
                  EXPOSED( i ) = j;
                }
              else if ( i != MATE( j ) )
                {
                  add_arc( &A, i, MATE( j ) );
                }
            }
        }
    }
}

/*
 * See Figure 11-2, "The Hungarian method", page 251.
 *
 * Corresponds to lines 26--27, 38--39.
 * Called by hm_pre_search() and hm_search().
 */
static void
hm_update_slack( hm_data *hm, int z )
{
  int k, tmp;
  for EACH_U( k )
    {
      tmp = C( z, k ) - ALPHA( z ) - BETA( k );
      if ( 0 <= tmp && tmp < SLACK( k ) )
        {
          SLACK( k ) = tmp;
          /*
           * The following decrement and increment are necessary to maintain
           * the count[] array, which is not included in the original Figure
           * 11-2 implementation, and whose addition and purpose are described
           * above in hm_construct_auxiliary_graph().
           */
          if ( NHBOR( k ) != blank )
            {
              --COUNT( NHBOR( k ) );
            }
          ++COUNT( z );
          NHBOR( k ) = z;
        }
    }
}

/*
 * See Figure 11-2, "The Hungarian method", page 251.
 *
 * Corresponds to lines 19--28.
 */
static bool
hm_pre_search( hm_data *hm )
{
  int i;
  Q.size = 0;
  for EACH_V( i )
    {
      if ( MATE( i ) == blank )
        {
          if ( EXPOSED( i ) != blank )
            {
              hm_augment( hm, i );
              return false; /* goto endstage */
            }
          stack_push( &Q, i );
          LABEL( i ) = blank;
          hm_update_slack( hm, i );
        }
    }
  return true;
}

/*
 * See Figure 11-2, "The Hungarian method", page 251.
 *
 * Corresponds to lines 29--41.
 */
static bool
hm_search( hm_data *hm )
{
  int i, j, z;
  while ( Q.size != 0 )
    {
      i = stack_pop( &Q );
      for ( z = 0; z < A.size; ++z )
        {
          if ( A.data[ z ].x == i )
            {
              j = A.data[ z ].y;
              if ( LABEL( j ) == blank )
                {
                  LABEL( j ) = i;
                  if ( EXPOSED( j ) != blank )
                    {
                      hm_augment( hm, j );
                      return false; /* goto endstage */
                    }
                  /*
                   * The following instruction is listed just before the prior
                   * conditional in Figure 11-2. Here, it is relocated simply
                   * because its execution would serve no purpose if the prior
                   * conditional executes.
                   */
                  stack_push( &Q, j );
                  hm_update_slack( hm, j );
                }
            }
        }
    }
  return true;
}

/*
 * See Figure 11-2, "The Hungarian method", page 252.
 *
 * Corresponds to "procedure modify".
 */
static bool
hm_modify( hm_data *hm )
{
  int i, j, theta_one;
  /*
   * Determine theta_one.
   */
  theta_one = INT_MAX;
  for EACH_U( j )
    {
      if ( 0 < SLACK( j ) && SLACK( j ) < theta_one )
        {
          theta_one = SLACK( j );
        }
    }
  theta_one /= 2;
  /*
   * Update the dual variable alpha.
   */
  for EACH_V( i )
    {
      /*
       * The following conditional expression has been changed from its form
       * in Figure 11-2. Here, an additional check on the count[] array is
       * performed to account for a certain type of "labeling" that is
       * mentioned in the Example 11.1 walk-through but is omitted from the
       * Figure 11-2 implementation.
       *
       * See the comments provided near the initialization of count[] in the
       * function hm_construct_auxiliary_graph().
       */
      if ( LABEL( i ) != blank || COUNT( i ) > 0 )
        {
          ALPHA( i ) += theta_one;
        }
      else
        {
          ALPHA( i ) -= theta_one;
        }
    }
  /*
   * Update the dual variable beta.
   */
  for EACH_U( j )
    {
      if ( SLACK( j ) == 0 )
        {
          BETA( j ) -= theta_one;
        }
      else
        {
          BETA( j ) += theta_one;
        }
    }
  /*
   * Update slack and check for new admissible edges.
   */
  for EACH_U( j )
    {
      if ( SLACK( j ) > 0 )
        {
          SLACK( j ) -= 2 * theta_one;
          if ( SLACK( j ) == 0 )
            {
              if ( MATE( j ) == blank )
                {
                  EXPOSED( NHBOR( j ) ) = j;
                  hm_augment( hm, NHBOR( j ) );
                  return false; /* goto endstage */
                }
              else
                {
                  /*
                   * The following statement corresponds to a pseudo-code
                   * command that should be removed from the else-clause of
                   * the modify procedure in Figure 11-2.
                   *
                   * LABEL( MATE( j ) ) = NHBOR( j );
                   *
                   * The inclusion of the above statement causes the arc
                   * added in one of the next statements to never be considered
                   * in following "search" sub-stages during this stage, and it
                   * partially duplicates what would happen in these sub-stages
                   * if the arc were to be considered there. The result of
                   * inclusion is (often) non-optimality of the algorithm's
                   * output.
                   */
                  /*
                   * The next statement corresponds to a pseudo-code command
                   * (in the same else-clause) that should be modified
                   * slightly. In Figure 11-2, this command "pushes" mate[ u ]
                   * into Q when it should be "pushing" nhbor[ u ] instead.
                   * This is because the purpose of this command is to ensure
                   * that the soon-to-be-added arc will be considered in the
                   * next "search" sub-stage, and consideration is dependent
                   * upon the arc-tail, not the arc-head.
                   */
                  stack_push( &Q, NHBOR( j ) ); /* Note modification */
                  add_arc( &A, NHBOR( j ), MATE( j ) );
                }
            }
        }
    }
  return true;
}

/*
 * See Figure 11-2, "The Hungarian method", pages 251--252.
 *
 * Input:
 *
 * mate  Points to a memory block of at least 2 * n ints. It is used to
 *       represent and return the solution matching. See page 223 for a
 *       contextual description.
 *
 * c     Points to a memory block of at least n * n ints. It contains the n*n
 *       cost matrix c[ 0...n-1 ][ 0...n-1 ] that implicitly defines the
 *       complete bipartite graph G=(V,U,E). The left and right indices
 *       respectively comprise vertex labels from V and U.
 *
 * n     Is the size of V and the size of U.
 *
 * Output:
 *
 * Returns mate filled with the correct values to represent the solution
 * matching, where V={0,...,n-1} and U={n,...,2n-1}. An edge (v,u) is part of
 * the matching if and only if (v,mate[v])=(mate[u],u).
 */
int *
hungarian_method( int *mate, int *c, int n )
{
  /*
   * Initialize the algorithm's internal data structures.
   */
  hm_data hm;
  if ( !hm_data_internal_allocate( &hm, n ) )
    return NULL;
  hm.mate = mate;
  hm.c = c;
  /*
   * Double each cost to ensure integrality of the alphabeta algorithm.
   */
  int i, j;
  for ( i = 0; i < n; ++i )
    for ( j = 0; j < n; ++j )
      c[ i * n + j ] *= 2;
  /*
   * Run the Hungarian method as described in Section 11.2 and Figure 11-2.
   */
  int s;
  hm_initialize( &hm );
  hm.q.size = 0;
  hm.a.size = 0;
  for ( s = 1; s <= n; ++s )
    {
      hm_construct_auxiliary_graph( &hm );
      if ( hm_pre_search( &hm ) )
        while( hm_search( &hm ) && hm_modify( &hm ) );
    }
  /*
   * Reset (halve) each cost and clean up.
   */
  for ( i = 0; i < n; ++i )
    for ( j = 0; j < n; ++j )
      c[ i * n + j ] /= 2;
  hm_data_internal_free( &hm );
  return mate;
}
