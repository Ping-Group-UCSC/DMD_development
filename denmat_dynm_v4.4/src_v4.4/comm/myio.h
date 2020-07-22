#pragma once
#include <sstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <string>
#include <scalar.h>
using namespace std;

void printf_complex_mat(complex *m, int n, string s);

void fprintf_complex_mat(FILE *fp, complex *m, int n, string s);

void printf_complex_mat(complex *a, int m, int n, string s);

void fprintf_complex_mat(FILE *fp, complex *a, int m, int n, string s);

void error_message(string s, string routine="");

void check_file_size(FILE *fp, size_t expect_size, string message);

void fseek_bigfile(FILE *fp, size_t count, size_t size);

inline bool exists(string name);

int last_file_index(string pre, string suf);

bool is_dir(string name);