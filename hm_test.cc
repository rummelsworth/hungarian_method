// Copyright 2010, 2017 William Rummler (w.a.rummler@gmail.com)
//
// This program is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with
// this program. If not, see <http://www.gnu.org/licenses/>.

// This file contains a very basic hard-coded harness for testing the
// implementation found in hungarian_method.c against the "brute force"
// implementation found in brute_force_assignment.c for solving the assignment
// problem.

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ctime>

#include "brute_force_assignment.h"
#include "hungarian_method.h"

#define TEST_DIM 8
#define NUM_TESTS 1000
#define MAX_COST 100

// This function will fill the n*n cost matrix c with random values between one
// and MAX_COST.
static void fill_randomly(int *c, int n)
{
    int i, j;
    for (i = 0; i < n; ++i)
    {
        for (j = 0; j < n; ++j)
        {
            c[i * n + j] = rand() % MAX_COST + 1;
        }
    }
}

// This function computes the cost of the complete matching represented by mate
// with respect to the n*n cost matrix c.
static int compute_cost(int *mate, int *c, int n)
{
    int i, cost = 0;
    for (i = 0; i < n; ++i)
    {
        cost += c[i * n + mate[i] - n];
    }

    return cost;
}

// This is a convenience function for displaying the contents of the n*n cost
// matrix c when tracing execution.
void print_c(int *c, int n)
{
    int i, j;
    for (i = 0; i < n; ++i)
    {
        for (j = 0; j < n; ++j)
        {
            printf(" %2d", c[i * n + j]);
        }

        printf("\n");
    }
}

// The test harness entry point.
int main(void)
{
    int test;

    // Initialized to [Papadimitriou & Steiglitz] Example 11.1's input, for basic
    // non-random testing when TEST_DIM is set to 5.
    int c[TEST_DIM * TEST_DIM] = { 7, 2, 1, 9, 4, \
                                   9, 6, 9, 5, 5, \
                                   3, 8, 3, 1, 8, \
                                   7, 9, 4, 2, 2, \
                                   8, 4, 7, 4, 8 };

    int mate[2 * TEST_DIM];
    int bf_cost, hm_cost;
    int num_pass = 0;
    srand((unsigned)time(NULL));
    for (test = 1; test <= NUM_TESTS; ++test)
    {
        // Prepare the test. Comment this out if testing the Example 11.1 input
        // to which c is initialized above. You may also want to change NUM_TESTS
        // to 1 in that case.
        fill_randomly(c, TEST_DIM);

        // Compute the best cost via brute-force.
        brute_force_assignment(mate, c, TEST_DIM);
        bf_cost = compute_cost(mate, c, TEST_DIM);

        // Compute the cost via Hungarian method.
        hungarian_method(mate, c, TEST_DIM);
        hm_cost = compute_cost(mate, c, TEST_DIM);

        // Check and display output.
        if (bf_cost == hm_cost)
        {
            ++num_pass;
        }

        printf(
            "Test %d: %s\n" \
            "         Brute Force = %10d\n" \
            "    Hungarian Method = %10d\n",
            test,
            bf_cost == hm_cost ? "+++ Pass +++" : "--- Fail ---",
            bf_cost, hm_cost);
    }

    printf("Number of tests passed = %d out of %d\n", num_pass, NUM_TESTS);
    return 0;
}
