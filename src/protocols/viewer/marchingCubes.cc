// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author

#include <protocols/viewer/marchingCubes.hh>

namespace protocols {
namespace viewer {


const int POLY_CASES[][21] = {
{ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0},
{ 1,  9,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 1 */
{ 1,  2, 10,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 2 */
{ 2,  9,  4,  0, 10,  9,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 3 */
{ 2,  3, 12,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 4 */
{ 1,  9,  4,  0,  2,  3, 12,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 5 */
{10,  3, 12,  0,  1,  3, 10,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 6 */
{ 3,  9,  4,  0,  3, 12,  9,  0, 12, 10,  9,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 7 */
{ 4, 11,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 8 */
{ 1, 11,  3,  0,  9, 11,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 9 */
{ 2, 10,  1,  0,  3,  4, 11,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 10 */
{ 2, 11,  3,  0,  2, 10, 11,  0, 10,  9, 11,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 11 */
{ 4, 12,  2,  0, 11, 12,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 12 */
{ 1, 12,  2,  0,  1,  9, 12,  0,  9, 11, 12,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 13 */
{ 4, 10,  1,  0,  4, 11, 10,  0, 11, 12, 10,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 14 */
{10,  9, 12,  0, 12,  9, 11,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 15 */
{ 5,  8,  9,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 16 */
{ 5,  4,  1,  0,  8,  4,  5,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 17 */
{ 1,  2, 10,  0,  9,  5,  8,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 18 */
{ 5,  2, 10,  0,  5,  8,  2,  0,  8,  4,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 19 */
{ 2,  3, 12,  0,  9,  5,  8,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 20 */
{ 4,  5,  8,  0,  4,  1,  5,  0,  2,  3, 12,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 21 */
{10,  3, 12,  0, 10,  1,  3,  0,  9,  5,  8,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 22 */
{ 3, 12, 10,  0,  3, 10,  8,  0,  3,  8,  4,  0,  8, 10,  5,  0,  0,  0,  0,  0, 0}, /* 23 */
{ 9,  5,  8,  0,  4, 11,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 24 */
{11,  5,  8,  0, 11,  3,  5,  0,  3,  1,  5,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 25 */
{10,  1,  2,  0,  9,  5,  8,  0,  3,  4, 11,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 26 */
{ 5,  8, 11,  0, 10,  5, 11,  0, 10, 11,  3,  0, 10,  3,  2,  0,  0,  0,  0,  0, 0}, /* 27 */
{ 4, 12,  2,  0,  4, 11, 12,  0,  8,  9,  5,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 28 */
{ 2, 11, 12,  0,  2,  5, 11,  0,  2,  1,  5,  0,  8, 11,  5,  0,  0,  0,  0,  0, 0}, /* 29 */
{ 5,  8,  9,  0, 10,  1, 11,  0, 10, 11, 12,  0, 11,  1,  4,  0,  0,  0,  0,  0, 0}, /* 30 */
{ 5,  8, 11,  0,  5, 11, 10,  0, 10, 11, 12,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 31 */
{10,  6,  5,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 32 */
{10,  6,  5,  0,  1,  9,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 33 */
{ 1,  6,  5,  0,  2,  6,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 34 */
{ 9,  6,  5,  0,  9,  4,  6,  0,  4,  2,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 35 */
{ 2,  3, 12,  0, 10,  6,  5,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 36 */
{ 4,  1,  9,  0,  2,  3, 12,  0,  5, 10,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 37 */
{ 6,  3, 12,  0,  6,  5,  3,  0,  5,  1,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 38 */
{ 3, 12,  6,  0,  4,  3,  6,  0,  4,  6,  5,  0,  4,  5,  9,  0,  0,  0,  0,  0, 0}, /* 39 */
{10,  6,  5,  0,  3,  4, 11,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 40 */
{ 1, 11,  3,  0,  1,  9, 11,  0,  5, 10,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 41 */
{ 1,  6,  5,  0,  1,  2,  6,  0,  3,  4, 11,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 42 */
{ 3,  2,  6,  0,  3,  6,  9,  0,  3,  9, 11,  0,  5,  9,  6,  0,  0,  0,  0,  0, 0}, /* 43 */
{12,  4, 11,  0, 12,  2,  4,  0, 10,  6,  5,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 44 */
{ 5, 10,  6,  0,  1,  9,  2,  0,  9, 12,  2,  0,  9, 11, 12,  0,  0,  0,  0,  0, 0}, /* 45 */
{ 6,  5,  1,  0,  6,  1, 11,  0,  6, 11, 12,  0, 11,  1,  4,  0,  0,  0,  0,  0, 0}, /* 46 */
{ 6,  5,  9,  0,  6,  9, 12,  0, 12,  9, 11,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 47 */
{10,  8,  9,  0,  6,  8, 10,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 48 */
{10,  4,  1,  0, 10,  6,  4,  0,  6,  8,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 49 */
{ 1,  8,  9,  0,  1,  2,  8,  0,  2,  6,  8,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 50 */
{ 2,  6,  4,  0,  4,  6,  8,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 51 */
{10,  8,  9,  0, 10,  6,  8,  0, 12,  2,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 52 */
{12,  2,  3,  0, 10,  6,  1,  0,  6,  4,  1,  0,  6,  8,  4,  0,  0,  0,  0,  0, 0}, /* 53 */
{ 9,  1,  3,  0,  9,  3,  6,  0,  9,  6,  8,  0, 12,  6,  3,  0,  0,  0,  0,  0, 0}, /* 54 */
{ 3, 12,  6,  0,  3,  6,  4,  0,  4,  6,  8,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 55 */
{ 8, 10,  6,  0,  8,  9, 10,  0,  4, 11,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 56 */
{10,  6,  8,  0, 10,  8,  3,  0, 10,  3,  1,  0,  3,  8, 11,  0,  0,  0,  0,  0, 0}, /* 57 */
{ 3,  4, 11,  0,  1,  2,  9,  0,  2,  8,  9,  0,  2,  6,  8,  0,  0,  0,  0,  0, 0}, /* 58 */
{11,  3,  2,  0, 11,  2,  8,  0,  8,  2,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 59 */
{10,  6,  9,  0,  9,  6,  8,  0, 12,  2,  4,  0, 12,  4, 11,  0,  0,  0,  0,  0, 0}, /* 60 */
{ 6,  8,  1,  0,  6,  1, 10,  0,  8, 11,  1,  0,  2,  1, 12,  0, 11, 12,  1,  0, 0}, /* 61 */
{11, 12,  1,  0, 11,  1,  4,  0, 12,  6,  1,  0,  9,  1,  8,  0,  6,  8,  1,  0, 0}, /* 62 */
{11, 12,  6,  0,  8, 11,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 63 */
{12,  7,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 64 */
{ 1,  9,  4,  0,  6, 12,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 65 */
{10,  1,  2,  0,  6, 12,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 66 */
{ 2,  9,  4,  0,  2, 10,  9,  0,  6, 12,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 67 */
{ 2,  7,  6,  0,  3,  7,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 68 */
{ 2,  7,  6,  0,  2,  3,  7,  0,  4,  1,  9,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 69 */
{10,  7,  6,  0, 10,  1,  7,  0,  1,  3,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 70 */
{ 6, 10,  9,  0,  6,  9,  3,  0,  6,  3,  7,  0,  4,  3,  9,  0,  0,  0,  0,  0, 0}, /* 71 */
{ 3,  4, 11,  0, 12,  7,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 72 */
{11,  1,  9,  0, 11,  3,  1,  0, 12,  7,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 73 */
{ 1,  2, 10,  0,  3,  4, 11,  0,  6, 12,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 74 */
{ 6, 12,  7,  0,  2, 10,  3,  0, 10, 11,  3,  0, 10,  9, 11,  0,  0,  0,  0,  0, 0}, /* 75 */
{ 7,  4, 11,  0,  7,  6,  4,  0,  6,  2,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 76 */
{ 1,  9, 11,  0,  1, 11,  6,  0,  1,  6,  2,  0,  6, 11,  7,  0,  0,  0,  0,  0, 0}, /* 77 */
{ 4, 11,  7,  0,  1,  4,  7,  0,  1,  7,  6,  0,  1,  6, 10,  0,  0,  0,  0,  0, 0}, /* 78 */
{ 7,  6, 10,  0,  7, 10, 11,  0, 11, 10,  9,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 79 */
{ 6, 12,  7,  0,  5,  8,  9,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 80 */
{ 5,  4,  1,  0,  5,  8,  4,  0,  7,  6, 12,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 81 */
{ 2, 10,  1,  0,  6, 12,  7,  0,  9,  5,  8,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 82 */
{12,  7,  6,  0,  2, 10,  8,  0,  2,  8,  4,  0,  8, 10,  5,  0,  0,  0,  0,  0, 0}, /* 83 */
{ 7,  2,  3,  0,  7,  6,  2,  0,  5,  8,  9,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 84 */
{ 2,  3,  6,  0,  6,  3,  7,  0,  4,  1,  5,  0,  4,  5,  8,  0,  0,  0,  0,  0, 0}, /* 85 */
{ 9,  5,  8,  0, 10,  1,  6,  0,  1,  7,  6,  0,  1,  3,  7,  0,  0,  0,  0,  0, 0}, /* 86 */
{ 8,  4, 10,  0,  8, 10,  5,  0,  4,  3, 10,  0,  6, 10,  7,  0,  3,  7, 10,  0, 0}, /* 87 */
{ 4, 11,  3,  0,  8,  9,  5,  0, 12,  7,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 88 */
{ 6, 12,  7,  0,  5,  8,  3,  0,  5,  3,  1,  0,  3,  8, 11,  0,  0,  0,  0,  0, 0}, /* 89 */
{ 1,  2, 10,  0,  5,  8,  9,  0,  3,  4, 11,  0,  6, 12,  7,  0,  0,  0,  0,  0, 0}, /* 90 */
{10,  3,  2,  0, 10, 11,  3,  0, 10,  5, 11,  0,  8, 11,  5,  0,  6, 12,  7,  0, 0}, /* 91 */
{ 9,  5,  8,  0,  4, 11,  6,  0,  4,  6,  2,  0,  6, 11,  7,  0,  0,  0,  0,  0, 0}, /* 92 */
{ 6,  2, 11,  0,  6, 11,  7,  0,  2,  1, 11,  0,  8, 11,  5,  0,  1,  5, 11,  0, 0}, /* 93 */
{ 1,  6, 10,  0,  1,  7,  6,  0,  1,  4,  7,  0, 11,  7,  4,  0,  9,  5,  8,  0, 0}, /* 94 */
{ 7,  6, 10,  0,  7, 10, 11,  0,  5,  8, 10,  0,  8, 11, 10,  0,  0,  0,  0,  0, 0}, /* 95 */
{12,  5, 10,  0,  7,  5, 12,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 96 */
{ 5, 12,  7,  0,  5, 10, 12,  0,  1,  9,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 97 */
{12,  1,  2,  0, 12,  7,  1,  0,  7,  5,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 98 */
{ 9,  4,  2,  0,  9,  2,  7,  0,  9,  7,  5,  0,  7,  2, 12,  0,  0,  0,  0,  0, 0}, /* 99  */
{ 2,  5, 10,  0,  2,  3,  5,  0,  3,  7,  5,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 100 */
{ 4,  1,  9,  0,  2,  3, 10,  0,  3,  5, 10,  0,  3,  7,  5,  0,  0,  0,  0,  0, 0}, /* 101 */
{ 1,  3,  5,  0,  5,  3,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 102 */
{ 9,  4,  3,  0,  9,  3,  5,  0,  5,  3,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 103 */
{12,  5, 10,  0, 12,  7,  5,  0, 11,  3,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 104 */
{ 1,  9,  3,  0,  3,  9, 11,  0,  5, 10, 12,  0,  5, 12,  7,  0,  0,  0,  0,  0, 0}, /* 105 */
{ 4, 11,  3,  0,  1,  2,  7,  0,  1,  7,  5,  0,  7,  2, 12,  0,  0,  0,  0,  0, 0}, /* 106 */
{ 7,  5,  2,  0,  7,  2, 12,  0,  5,  9,  2,  0,  3,  2, 11,  0,  9, 11,  2,  0, 0}, /* 107 */
{10,  7,  5,  0, 10,  4,  7,  0, 10,  2,  4,  0, 11,  7,  4,  0,  0,  0,  0,  0, 0}, /* 108 */
{ 9, 11,  2,  0,  9,  2,  1,  0, 11,  7,  2,  0, 10,  2,  5,  0,  7,  5,  2,  0, 0}, /* 109 */
{ 4, 11,  7,  0,  4,  7,  1,  0,  1,  7,  5,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 110 */
{ 7,  5,  9,  0, 11,  7,  9,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 111 */
{ 8, 12,  7,  0,  8,  9, 12,  0,  9, 10, 12,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 112 */
{ 1,  8,  4,  0,  1, 12,  8,  0,  1, 10, 12,  0,  7,  8, 12,  0,  0,  0,  0,  0, 0}, /* 113 */
{12,  7,  8,  0,  2, 12,  8,  0,  2,  8,  9,  0,  2,  9,  1,  0,  0,  0,  0,  0, 0}, /* 114 */
{12,  7,  8,  0, 12,  8,  2,  0,  2,  8,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 115 */
{ 2,  3,  7,  0,  2,  7,  9,  0,  2,  9, 10,  0,  9,  7,  8,  0,  0,  0,  0,  0, 0}, /* 116 */
{ 3,  7, 10,  0,  3, 10,  2,  0,  7,  8, 10,  0,  1, 10,  4,  0,  8,  4, 10,  0, 0}, /* 117 */
{ 8,  9,  1,  0,  8,  1,  7,  0,  7,  1,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 118 */
{ 8,  4,  3,  0,  7,  8,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 119 */
{ 3,  4, 11,  0, 12,  7,  9,  0, 12,  9, 10,  0,  9,  7,  8,  0,  0,  0,  0,  0, 0}, /* 120 */
{ 3,  1,  8,  0,  3,  8, 11,  0,  1, 10,  8,  0,  7,  8, 12,  0, 10, 12,  8,  0, 0}, /* 121 */
{ 2,  9,  1,  0,  2,  8,  9,  0,  2, 12,  8,  0,  7,  8, 12,  0,  3,  4, 11,  0, 0}, /* 122 */
{11,  3,  2,  0, 11,  2,  8,  0, 12,  7,  2,  0,  7,  8,  2,  0,  0,  0,  0,  0, 0}, /* 123 */
{ 9, 10,  7,  0,  9,  7,  8,  0, 10,  2,  7,  0, 11,  7,  4,  0,  2,  4,  7,  0, 0}, /* 124 */
{ 1, 10,  2,  0, 11,  7,  8,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 125 */
{ 8,  9,  1,  0,  8,  1,  7,  0,  4, 11,  1,  0, 11,  7,  1,  0,  0,  0,  0,  0, 0}, /* 126 */
{ 8, 11,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 127 */
{ 8,  7, 11,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 128 */
{ 4,  1,  9,  0, 11,  8,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 129 */
{ 1,  2, 10,  0, 11,  8,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 130 */
{ 9,  2, 10,  0,  9,  4,  2,  0, 11,  8,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 131 */
{12,  2,  3,  0,  7, 11,  8,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 132 */
{ 2,  3, 12,  0,  4,  1,  9,  0,  7, 11,  8,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 133 */
{ 3, 10,  1,  0,  3, 12, 10,  0,  7, 11,  8,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 134 */
{ 7, 11,  8,  0,  3, 12,  4,  0, 12,  9,  4,  0, 12, 10,  9,  0,  0,  0,  0,  0, 0}, /* 135 */
{ 8,  3,  4,  0,  7,  3,  8,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 136 */
{ 8,  1,  9,  0,  8,  7,  1,  0,  7,  3,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 137 */
{ 3,  8,  7,  0,  3,  4,  8,  0,  1,  2, 10,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 138 */
{ 2,  7,  3,  0,  2,  9,  7,  0,  2, 10,  9,  0,  9,  8,  7,  0,  0,  0,  0,  0, 0}, /* 139 */
{12,  8,  7,  0, 12,  2,  8,  0,  2,  4,  8,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 140 */
{12,  8,  7,  0,  2,  8, 12,  0,  2,  9,  8,  0,  2,  1,  9,  0,  0,  0,  0,  0, 0}, /* 141 */
{ 1,  4,  8,  0,  1,  8, 12,  0,  1, 12, 10,  0,  7, 12,  8,  0,  0,  0,  0,  0, 0}, /* 142 */
{ 8,  7, 12,  0,  8, 12,  9,  0,  9, 12, 10,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 143 */
{ 7,  9,  5,  0, 11,  9,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 144 */
{ 4,  7, 11,  0,  4,  1,  7,  0,  1,  5,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 145 */
{ 9,  7, 11,  0,  9,  5,  7,  0, 10,  1,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 146 */
{10,  5,  7,  0, 10,  7,  4,  0, 10,  4,  2,  0, 11,  4,  7,  0,  0,  0,  0,  0, 0}, /* 147 */
{ 7,  9,  5,  0,  7, 11,  9,  0,  3, 12,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 148 */
{ 2,  3, 12,  0,  4,  1, 11,  0,  1,  7, 11,  0,  1,  5,  7,  0,  0,  0,  0,  0, 0}, /* 149 */
{ 5, 11,  9,  0,  5,  7, 11,  0,  1,  3, 10,  0,  3, 12, 10,  0,  0,  0,  0,  0, 0}, /* 150 */
{12, 10,  4,  0, 12,  4,  3,  0, 10,  5,  4,  0, 11,  4,  7,  0,  5,  7,  4,  0, 0}, /* 151 */
{ 9,  3,  4,  0,  9,  5,  3,  0,  5,  7,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 152 */
{ 1,  5,  3,  0,  5,  7,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 153 */
{ 2, 10,  1,  0,  3,  4,  5,  0,  3,  5,  7,  0,  5,  4,  9,  0,  0,  0,  0,  0, 0}, /* 154 */
{ 2, 10,  5,  0,  2,  5,  3,  0,  3,  5,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 155 */
{ 9,  2,  4,  0,  9,  7,  2,  0,  9,  5,  7,  0,  7, 12,  2,  0,  0,  0,  0,  0, 0}, /* 156 */
{12,  2,  1,  0, 12,  1,  7,  0,  7,  1,  5,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 157 */
{ 5,  7,  4,  0,  5,  4,  9,  0,  7, 12,  4,  0,  1,  4, 10,  0, 12, 10,  4,  0, 0}, /* 158 */
{12, 10,  5,  0,  7, 12,  5,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 159 */
{ 5, 10,  6,  0,  8,  7, 11,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 160 */
{ 1,  9,  4,  0,  5, 10,  6,  0, 11,  8,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 161 */
{ 6,  1,  2,  0,  6,  5,  1,  0,  8,  7, 11,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 162 */
{11,  8,  7,  0,  9,  4,  5,  0,  4,  6,  5,  0,  4,  2,  6,  0,  0,  0,  0,  0, 0}, /* 163 */
{10,  6,  5,  0, 12,  2,  3,  0,  8,  7, 11,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 164 */
{ 7, 11,  8,  0,  2,  3, 12,  0,  1,  9,  4,  0,  5, 10,  6,  0,  0,  0,  0,  0, 0}, /* 165 */
{ 8,  7, 11,  0,  6,  5, 12,  0,  5,  3, 12,  0,  5,  1,  3,  0,  0,  0,  0,  0, 0}, /* 166 */
{ 4,  5,  9,  0,  4,  6,  5,  0,  4,  3,  6,  0, 12,  6,  3,  0, 11,  8,  7,  0, 0}, /* 167 */
{ 8,  3,  4,  0,  8,  7,  3,  0,  6,  5, 10,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 168 */
{10,  6,  5,  0,  1,  9,  7,  0,  1,  7,  3,  0,  7,  9,  8,  0,  0,  0,  0,  0, 0}, /* 169 */
{ 4,  7,  3,  0,  4,  8,  7,  0,  2,  6,  1,  0,  6,  5,  1,  0,  0,  0,  0,  0, 0}, /* 170 */
{ 7,  3,  9,  0,  7,  9,  8,  0,  3,  2,  9,  0,  5,  9,  6,  0,  2,  6,  9,  0, 0}, /* 171 */
{10,  6,  5,  0, 12,  2,  7,  0,  2,  8,  7,  0,  2,  4,  8,  0,  0,  0,  0,  0, 0}, /* 172 */
{ 2,  7, 12,  0,  2,  8,  7,  0,  2,  1,  8,  0,  9,  8,  1,  0, 10,  6,  5,  0, 0}, /* 173 */
{ 5,  1, 12,  0,  5, 12,  6,  0,  1,  4, 12,  0,  7, 12,  8,  0,  4,  8, 12,  0, 0}, /* 174 */
{ 8,  7, 12,  0,  8, 12,  9,  0,  6,  5, 12,  0,  5,  9, 12,  0,  0,  0,  0,  0, 0}, /* 175 */
{ 7, 10,  6,  0,  7, 11, 10,  0, 11,  9, 10,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 176 */
{ 4,  7, 11,  0,  1,  7,  4,  0,  1,  6,  7,  0,  1, 10,  6,  0,  0,  0,  0,  0, 0}, /* 177 */
{ 1, 11,  9,  0,  1,  6, 11,  0,  1,  2,  6,  0,  6,  7, 11,  0,  0,  0,  0,  0, 0}, /* 178 */
{ 7, 11,  4,  0,  7,  4,  6,  0,  6,  4,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 179 */
{ 2,  3, 12,  0, 10,  6, 11,  0, 10, 11,  9,  0, 11,  6,  7,  0,  0,  0,  0,  0, 0}, /* 180 */
{ 1, 11,  4,  0,  1,  7, 11,  0,  1, 10,  7,  0,  6,  7, 10,  0,  2,  3, 12,  0, 0}, /* 181 */
{11,  9,  6,  0, 11,  6,  7,  0,  9,  1,  6,  0, 12,  6,  3,  0,  1,  3,  6,  0, 0}, /* 182 */
{ 7, 11,  4,  0,  7,  4,  6,  0,  3, 12,  4,  0, 12,  6,  4,  0,  0,  0,  0,  0, 0}, /* 183 */
{ 6,  9, 10,  0,  6,  3,  9,  0,  6,  7,  3,  0,  4,  9,  3,  0,  0,  0,  0,  0, 0}, /* 184 */
{10,  6,  7,  0, 10,  7,  1,  0,  1,  7,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 185 */
{ 2,  6,  9,  0,  2,  9,  1,  0,  6,  7,  9,  0,  4,  9,  3,  0,  7,  3,  9,  0, 0}, /* 186 */
{ 2,  6,  7,  0,  3,  2,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 187 */
{ 2,  4,  7,  0,  2,  7, 12,  0,  4,  9,  7,  0,  6,  7, 10,  0,  9, 10,  7,  0, 0}, /* 188 */
{12,  2,  1,  0, 12,  1,  7,  0, 10,  6,  1,  0,  6,  7,  1,  0,  0,  0,  0,  0, 0}, /* 189 */
{ 1,  4,  9,  0,  6,  7, 12,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 190 */
{12,  6,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 191 */
{11,  6, 12,  0,  8,  6, 11,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 192 */
{11,  6, 12,  0, 11,  8,  6,  0,  9,  4,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 193 */
{ 6, 11,  8,  0,  6, 12, 11,  0,  2, 10,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 194 */
{12,  8,  6,  0, 12, 11,  8,  0, 10,  9,  2,  0,  9,  4,  2,  0,  0,  0,  0,  0, 0}, /* 195 */
{11,  2,  3,  0, 11,  8,  2,  0,  8,  6,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 196 */
{ 1,  9,  4,  0,  2,  3,  8,  0,  2,  8,  6,  0,  8,  3, 11,  0,  0,  0,  0,  0, 0}, /* 197 */
{10,  8,  6,  0, 10,  3,  8,  0, 10,  1,  3,  0,  3, 11,  8,  0,  0,  0,  0,  0, 0}, /* 198 */
{ 8,  6,  3,  0,  8,  3, 11,  0,  6, 10,  3,  0,  4,  3,  9,  0, 10,  9,  3,  0, 0}, /* 199 */
{ 3,  6, 12,  0,  3,  4,  6,  0,  4,  8,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 200 */
{ 9,  3,  1,  0,  9,  6,  3,  0,  9,  8,  6,  0, 12,  3,  6,  0,  0,  0,  0,  0, 0}, /* 201 */
{10,  1,  2,  0,  6, 12,  4,  0,  6,  4,  8,  0,  4, 12,  3,  0,  0,  0,  0,  0, 0}, /* 202 */
{10,  9,  3,  0, 10,  3,  2,  0,  9,  8,  3,  0, 12,  3,  6,  0,  8,  6,  3,  0, 0}, /* 203 */
{ 2,  4,  6,  0,  4,  8,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 204 */
{ 1,  9,  8,  0,  1,  8,  2,  0,  2,  8,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 205 */
{10,  1,  4,  0, 10,  4,  6,  0,  6,  4,  8,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 206 */
{10,  9,  8,  0,  6, 10,  8,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 207 */
{ 6,  9,  5,  0,  6, 12,  9,  0, 12, 11,  9,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 208 */
{ 6,  1,  5,  0,  6, 11,  1,  0,  6, 12, 11,  0, 11,  4,  1,  0,  0,  0,  0,  0, 0}, /* 209 */
{ 1,  2, 10,  0,  9,  5, 12,  0,  9, 12, 11,  0, 12,  5,  6,  0,  0,  0,  0,  0, 0}, /* 210 */
{12, 11,  5,  0, 12,  5,  6,  0, 11,  4,  5,  0, 10,  5,  2,  0,  4,  2,  5,  0, 0}, /* 211 */
{ 3,  6,  2,  0,  3,  9,  6,  0,  3, 11,  9,  0,  5,  6,  9,  0,  0,  0,  0,  0, 0}, /* 212 */
{ 1,  5, 11,  0,  1, 11,  4,  0,  5,  6, 11,  0,  3, 11,  2,  0,  6,  2, 11,  0, 0}, /* 213 */
{ 1,  3,  6,  0,  1,  6, 10,  0,  3, 11,  6,  0,  5,  6,  9,  0, 11,  9,  6,  0, 0}, /* 214 */
{10,  5,  6,  0,  3, 11,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 215 */
{ 3,  6, 12,  0,  4,  6,  3,  0,  4,  5,  6,  0,  4,  9,  5,  0,  0,  0,  0,  0, 0}, /* 216 */
{ 6, 12,  3,  0,  6,  3,  5,  0,  5,  3,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 217 */
{ 4, 12,  3,  0,  4,  6, 12,  0,  4,  9,  6,  0,  5,  6,  9,  0,  1,  2, 10,  0, 0}, /* 218 */
{ 6, 12,  3,  0,  6,  3,  5,  0,  2, 10,  3,  0, 10,  5,  3,  0,  0,  0,  0,  0, 0}, /* 219 */
{ 9,  5,  6,  0,  9,  6,  4,  0,  4,  6,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 220 */
{ 1,  5,  6,  0,  2,  1,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 221 */
{ 9,  5,  6,  0,  9,  6,  4,  0, 10,  1,  6,  0,  1,  4,  6,  0,  0,  0,  0,  0, 0}, /* 222 */
{10,  5,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 223 */
{ 5, 11,  8,  0,  5, 10, 11,  0, 10, 12, 11,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 224 */
{ 1,  9,  4,  0,  5, 10,  8,  0, 10, 11,  8,  0, 10, 12, 11,  0,  0,  0,  0,  0, 0}, /* 225 */
{ 2, 12, 11,  0,  2, 11,  5,  0,  2,  5,  1,  0,  8,  5, 11,  0,  0,  0,  0,  0, 0}, /* 226 */
{ 4,  2,  5,  0,  4,  5,  9,  0,  2, 12,  5,  0,  8,  5, 11,  0, 12, 11,  5,  0, 0}, /* 227 */
{ 5, 11,  8,  0, 10, 11,  5,  0, 10,  3, 11,  0, 10,  2,  3,  0,  0,  0,  0,  0, 0}, /* 228 */
{10,  8,  5,  0, 10, 11,  8,  0, 10,  2, 11,  0,  3, 11,  2,  0,  1,  9,  4,  0, 0}, /* 229 */
{11,  8,  5,  0, 11,  5,  3,  0,  3,  5,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 230 */
{11,  8,  5,  0, 11,  5,  3,  0,  9,  4,  5,  0,  4,  3,  5,  0,  0,  0,  0,  0, 0}, /* 231 */
{ 3, 10, 12,  0,  3,  8, 10,  0,  3,  4,  8,  0,  8,  5, 10,  0,  0,  0,  0,  0, 0}, /* 232 */
{10, 12,  8,  0, 10,  8,  5,  0, 12,  3,  8,  0,  9,  8,  1,  0,  3,  1,  8,  0, 0}, /* 233 */
{ 4,  8, 12,  0,  4, 12,  3,  0,  8,  5, 12,  0,  2, 12,  1,  0,  5,  1, 12,  0, 0}, /* 234 */
{ 2, 12,  3,  0,  9,  8,  5,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 235 */
{ 5, 10,  2,  0,  5,  2,  8,  0,  8,  2,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 236 */
{ 5, 10,  2,  0,  5,  2,  8,  0,  1,  9,  2,  0,  9,  8,  2,  0,  0,  0,  0,  0, 0}, /* 237 */
{ 5,  1,  4,  0,  8,  5,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 238 */
{ 5,  9,  8,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 239 */
{10, 12,  9,  0, 12, 11,  9,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 240 */
{ 4,  1, 10,  0,  4, 10, 11,  0, 11, 10, 12,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 241 */
{ 1,  2, 12,  0,  1, 12,  9,  0,  9, 12, 11,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 242 */
{ 4,  2, 12,  0, 11,  4, 12,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 243 */
{ 2,  3, 11,  0,  2, 11, 10,  0, 10, 11,  9,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 244 */
{ 4,  1, 10,  0,  4, 10, 11,  0,  2,  3, 10,  0,  3, 11, 10,  0,  0,  0,  0,  0, 0}, /* 245 */
{ 1,  3, 11,  0,  9,  1, 11,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 246 */
{ 4,  3, 11,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 247 */
{ 3,  4,  9,  0,  3,  9, 12,  0, 12,  9, 10,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 248 */
{10, 12,  3,  0,  1, 10,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 249 */
{ 3,  4,  9,  0,  3,  9, 12,  0,  1,  2,  9,  0,  2, 12,  9,  0,  0,  0,  0,  0, 0}, /* 250 */
{ 2, 12,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 251 */
{ 2,  4,  9,  0, 10,  2,  9,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 252 */
{ 1, 10,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 253 */
{ 1,  4,  9,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}, /* 254 */
{ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0}
};

const int VERTEX_OFF[][3] = {
{-1, -1, -1}, // dummy, no vertex with index 0
{0, 0, 0},
{1, 0, 0},
{1, 1, 0},
{0, 1, 0},
{0, 0, 1},
{1, 0, 1},
{1, 1, 1},
{0, 1, 1}
};

// Column labeling:
//
// Neighbor Vertex 1, Neighbor Vertex2, Direction.
//
// Direction is the vector from Neighbor 1 to Neighbor 2:
//   0 - +x
//   1 - +x
//   2 - +x
//   3 - -x
//   4 - -y
//   5 - -z
const int EDGE_NGHBRS[][3] = {
{-1, -1, -1}, // dummy, no edge with index 0

// v1  v2   d
{ 1,  2,  0}, //  1
{ 2,  3,  1}, //  2
{ 3,  4,  3}, //  3
{ 4,  1,  4}, //  4

{ 5,  6,  0}, //  5
{ 6,  7,  1}, //  6
{ 7,  8,  3}, //  7
{ 8,  5,  4}, //  8

{ 1,  5,  2}, //  9
{ 2,  6,  2}, // 10
{ 4,  8,  2}, // 11
{ 3,  7,  2}  // 12
};

}
}
