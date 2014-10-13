// ===========================================================
//     _/_/_/   _/_/_/  _/_/_/_/    _/_/_/_/  _/_/_/   _/_/_/
//      _/    _/       _/             _/    _/    _/   _/   _/
//     _/    _/       _/_/_/_/       _/    _/    _/   _/_/_/
//    _/    _/       _/             _/    _/    _/   _/
// _/_/_/   _/_/_/  _/_/_/_/_/     _/     _/_/_/   _/_/
// ===========================================================
//
// Common.h: the header file of SeqArray C/C++ codes
//
// Copyright (C) 2013 - 2014	Xiuwen Zheng [zhengx@u.washington.edu]
//
// This file is part of SeqArray.
//
// SeqArray is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License Version 3 as
// published by the Free Software Foundation.
//
// SeqArray is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with SeqArray.
// If not, see <http://www.gnu.org/licenses/>.


#include <R_GDS_CPP.h>
#include <dTrait.h>

#include <string>
#include <vector>
#include <list>
#include <map>

#include <cctype>
#include <cstring>


using namespace std;
using namespace CoreArray;


#define LongBool int



// ###########################################################
// The initialized object
// ###########################################################

/// the initial data
class COREARRAY_DLL_LOCAL TInitObject
{
public:
	struct TSelection
	{
		vector<C_BOOL> Sample;
		vector<C_BOOL> Variant;
	};

	typedef list<TSelection> TSelList;

	TInitObject();

	TSelection &Selection(SEXP gds);
	void Check_TrueArray(int Cnt);

	vector<C_BOOL> TRUE_ARRAY;
	vector<C_UInt8> GENO_BUFFER;
	map<int, TSelList> _Map;
};

extern TInitObject Init;



// ###########################################################
// define exception
// ###########################################################

class ErrSeqArray: public ErrCoreArray
{
public:
	ErrSeqArray(): ErrCoreArray()
		{ }
	ErrSeqArray(const char *fmt, ...): ErrCoreArray()
		{ _COREARRAY_ERRMACRO_(fmt); }
	ErrSeqArray(const std::string &msg): ErrCoreArray()
		{ fMessage = msg; }
};



// ###########################################################
// private functions
// ###########################################################

/// get the list element named str, or return NULL
COREARRAY_INLINE static SEXP getListElement(SEXP list, const char *str)
{
	SEXP elmt = R_NilValue;
	SEXP names = getAttrib(list, R_NamesSymbol);
	for (R_len_t i = 0; i < XLENGTH(list); i++)
	{
		if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0)
		{
			elmt = VECTOR_ELT(list, i);
			break;
		}
	}
	return elmt;
}

/// check CoreArray function
COREARRAY_INLINE static const char *SKIP(const char *p)
{
	while (isspace(*p)) p ++;
	return p;
}

/// check CoreArray function
COREARRAY_INLINE static string SHORT_TEXT(const char *p, int MaxNum=16)
{
	if ((int)strlen(p) <= MaxNum)
		return string(p);
	else
		return string(p, MaxNum) + "...";
}

/// get PdGDSObj from a SEXP object
COREARRAY_INLINE static void GDS_PATH_PREFIX_CHECK(const char *path)
{
	for (; *path != 0; path++)
	{
		if ((*path == '~') || (*path == '@'))
		{
			throw ErrSeqArray(
				"the variable name contains an invalid prefix '%c'.",
				*path);
		}
	}
}

/// check variable name
COREARRAY_INLINE static void GDS_VARIABLE_NAME_CHECK(const char *p)
{
	for (; *p != 0; p++)
	{
		if ((*p == '~') || (*p == '@') || (*p == '/'))
		{
			throw ErrSeqArray(
				"the variable name contains an invalid prefix '%c'.", *p);
		}
	}
}

/// get PdGDSObj from a SEXP object
COREARRAY_INLINE static string GDS_PATH_PREFIX(const string &path, char prefix)
{
	string s = path;
	for (int i=s.size()-1; i >= 0; i--)
	{
		if (s[i] == '/')
		{
			s.insert(i+1, &prefix, 1);
			return s;
		}
	}
	s.insert(s.begin(), prefix);
	return s;
}

/// get PdGDSObj from a SEXP object
COREARRAY_INLINE static string GDS_UP_PATH(const char *path)
{
	const char *p = path + strlen(path) - 1;
	while ((p!=path) && (*p != '/')) p --;
	return string(path, p);
}
