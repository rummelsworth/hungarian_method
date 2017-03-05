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

#include "brute_force_assignment.h"
#include <stdlib.h>
#include <string.h>
#include <limits.h>

/*
 * rho[] is an array used by the function perm_lex_successor() found below.
 * It is declared and maintained globally in order to be efficiently re-used
 * by numerous successive calls to that function.
 */
static int *rho;

/*
 * This function is registered with atexit() in order to free any memory
 * allocated for rho[] at program's end.
 */
static void
rho_free( void )
{
  if ( rho ) free( rho );
}

/*
 * An implementation of:
 *
 * Kreher & Stinson
 * "Combinatorial Algorithms: Generation, Enumeration, and Search"
 * Algorithm 2.14
 *
 * Generates the lexicographic successor of the given permutation.
 * Some lines are adjusted for zero-based indexing. Such lines are indicated.
 */
static int *
perm_lex_successor( int n, int *pi )
{
  static int m;
  int h, i, j, t, *p;
  i = n - 2; /* adjusted */
  while ( i >= 0 && pi[ i ] > pi[ i + 1 ] ) /* adjusted */
    {
      --i;
    }
  if ( i < 0 ) /* adjusted */
    return NULL;
  /*
   * BEGIN: memory management for rho[] (see declaration above)
   */
  if ( m < n )
    {
      if ( m == 0 )
        {
          atexit( rho_free );
        }
      p = realloc( rho, ( m = n ) * sizeof( *p ) );
      if ( !p )
        {
          return NULL;
        }
      rho = p;
    }
  /*
   * END: memory management for rho[]
   */
  j = n - 1; /* adjusted */
  while ( pi[ j ] < pi[ i ] )
    {
      --j;
    }
  t = pi[ j ];
  pi[ j ] = pi[ i ];
  pi[ i ] = t;
  for ( h = i + 1; h < n; ++h ) /* adjusted */
    {
      rho[ h ] = pi[ h ];
    }
  for ( h = i + 1; h < n; ++h ) /* adjusted */
    {
      pi[ h ] = rho[ n + i - h ]; /* adjusted */
    }
  return pi;
}

/*
 * This function uses an exhaustive "brute force" search over all possible
 * matchings to solve the assignment problem defined by the weighted complete
 * bipartite graph G=(V,U,E) which, in turn, is implicitly defined by the cost
 * matrix c given as an argument to this function.
 *
 * Specifically, it searches all n! permutations of {0,...,n-1} where
 * n=|V|=|U| and \pi_i is the i-th element of the permutation \pi, and
 * vertex i from V is assigned to vertex \pi_i from U.
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
brute_force_assignment( int *mate, int *c, int n )
{
  /*
   * Allocate temporary storage for exhaustive search.
   */
  int best_cost, current_cost;
  int *best_perm, *current_perm;
  if ( !( best_perm = malloc( n * sizeof( *best_perm ) ) )
       || !( current_perm = malloc( n * sizeof( *current_perm ) ) ) )
    return NULL;
  /*
   * Initialize the current permutation data structure with the
   * lexicographically first permutation [0,...,n-1].
   */
  int i;
  for ( i = 0; i < n; ++i )
    current_perm[ i ] = i;
  /*
   * Search over all n! permutations.
   */
  best_cost = INT_MAX;
  do
    {
      current_cost = 0;
      for ( i = 0; i < n; ++i )
        current_cost += c[ i * n + current_perm[ i ] ];
      if ( current_cost < best_cost )
        {
          best_cost = current_cost;
          memcpy( best_perm, current_perm, n * sizeof( *current_perm ) );
        }
    }
  while ( perm_lex_successor( n, current_perm ) );
  /*
   * Translate the best permutation into mate[].
   */
  int v, u;
  for ( v = 0; v < n; ++v )
    {
      u = n + best_perm[ v ];
      mate[ v ] = u;
      mate[ u ] = v;
    }
  return mate;
}
