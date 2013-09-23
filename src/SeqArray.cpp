// ===========================================================
//     _/_/_/   _/_/_/  _/_/_/_/    _/_/_/_/  _/_/_/   _/_/_/
//      _/    _/       _/             _/    _/    _/   _/   _/
//     _/    _/       _/_/_/_/       _/    _/    _/   _/_/_/
//    _/    _/       _/             _/    _/    _/   _/
// _/_/_/   _/_/_/  _/_/_/_/_/     _/     _/_/_/   _/_/
// ===========================================================
//
// main.cpp: SeqArray
//
// Copyright (C) 2013	Xiuwen Zheng
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

#include <CoreGDSLink.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <set>
#include <algorithm>

#include <R.h>
#include <Rdefines.h>
#include <stdio.h>

// R_XLEN_T_MAX is defined, R >= v3.0
#ifndef R_XLEN_T_MAX
#  define R_xlen_t	R_len_t
#  define XLENGTH	Rf_length
#endif


using namespace std;
using namespace CoreArray;
using namespace GDSInterface;



#define LongBool int


#ifdef COREARRAY_GNUG
#  ifdef COREARRAY_WINDOWS
#    define DLLEXPORT __attribute__((dllexport))
#  else
#    define DLLEXPORT
#  endif
#else
#  define DLLEXPORT __declspec(dllexport)
#endif



// ###########################################################
// defined macro
// ###########################################################

#define CORETRY			try {
#define CORETRY_CALL	bool has_error = false; CORETRY

#define CORECATCH(cmd)	} \
	catch (exception &E) { \
		gds_LastError() = E.what(); \
		cmd; \
	} \
	catch (const char *E) { \
		gds_LastError() = E; \
		cmd; \
	} \
	catch (...) { \
		gds_LastError() = "unknown error!"; \
		cmd; \
	}
#define CORECATCH_CALL	CORECATCH(has_error = true); \
	if (has_error) error(gds_LastError().c_str());





// ###########################################################
// define 
// ###########################################################

static const string BlackString;



// ###########################################################
// define exception
// ###########################################################

class ErrSeqArray: public ErrCoreArray
{
public:
	ErrSeqArray() {};
	ErrSeqArray(const char *fmt, ...) { _COREARRAY_ERRMACRO_(fmt); }
	ErrSeqArray(const std::string &msg) { fMessage = msg; }
};



// ###########################################################
// private functions
// ###########################################################

/// check CoreArray function
COREARRAY_INLINE static void CHECK(bool retval)
{
	if (!retval)
		throw ErrSeqArray(gds_LastError());
}
/// check CoreArray function
COREARRAY_INLINE static void* CHECK(void* retval)
{
	if (retval == NULL)
		throw ErrSeqArray(gds_LastError());
	return retval;
}
/// check CoreArray function
COREARRAY_INLINE static Int64 CHECK(Int64 retval)
{
	if (retval < 0)
		throw ErrSeqArray(gds_LastError());
	return retval;
}
/// check CoreArray function
COREARRAY_INLINE static size_t CHECK_SIZE(size_t retval)
{
	if (retval == ((size_t)-1))
		throw ErrSeqArray(gds_LastError());
	return retval;
}


/// check CoreArray function
COREARRAY_INLINE const char *SKIP(const char *p)
{
	while (isspace(*p)) p ++;
	return p;
}

/// check CoreArray function
COREARRAY_INLINE string SHORT_TEXT(const char *p, int MaxNum=16)
{
	if ((int)strlen(p) <= MaxNum)
		return string(p);
	else
		return string(p, MaxNum) + "...";
}

/// get PdGDSObj from a SEXP object
COREARRAY_INLINE PdGDSObj GDS_OBJECT(SEXP obj)
{
	PdGDSObj N;
	memcpy(&N, INTEGER(obj), sizeof(N));
	return N;
}

/// get PdGDSObj from a SEXP object
void GDS_PATH_PREFIX_CHECK(const char *path)
{
	for (; *path != 0; path++)
	{
		if ((*path == '~') || (*path == '@'))
			throw ErrSeqArray("the variable name contains an invalid prefix '%c'.", *path);
	}
}

/// check variable name
void GDS_VARIABLE_NAME_CHECK(const char *p)
{
	for (; *p != 0; p++)
	{
		if ((*p == '~') || (*p == '@') || (*p == '/'))
			throw ErrSeqArray("the variable name contains an invalid prefix '%c'.", *p);
	}
}

/// get PdGDSObj from a SEXP object
string GDS_PATH_PREFIX(const string &path, char prefix)
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
string GDS_UP_PATH(const char *path)
{
	const char *p = path + strlen(path) - 1;
	while ((p!=path) && (*p != '/')) p --;
	return string(path, p);
}

/// convert _SEXP to SEXP
COREARRAY_INLINE static SEXP _(_SEXP_ v)
{
	union {
		_SEXP_ f;
		SEXP t;
	} u;
	u.f = v;
	return u.t;
}

/// get the list element named str, or return NULL
static SEXP getListElement(SEXP list, const char *str)
{
	SEXP elmt = R_NilValue;
	SEXP names = getAttrib(list, R_NamesSymbol);
	for (R_len_t i = 0; i < length(list); i++)
	{
		if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0)
		{
			elmt = VECTOR_ELT(list, i);
			break;
		}
	}
	return elmt;
}



// ###########################################################
// The initialized object
// ###########################################################

/// the initial data
class TInitObject
{
public:
	struct TSelection
	{
		vector<CBOOL> Sample;
		vector<CBOOL> Variant;
	};

	typedef list<TSelection> TSelList;

	TInitObject(): TRUE_ARRAY(256, TRUE), GENO_BUFFER(1024) {}

	COREARRAY_INLINE TSelection &Selection(SEXP gds)
	{
		// TODO: check whether handle is valid
		int id = INTEGER(getListElement(gds, "id"))[0];
		TSelList &m = _Map[id];
		if (m.empty()) m.push_back(TSelection());
		return m.back();
	}

	COREARRAY_INLINE void Check_TrueArray(int Cnt)
	{
		if (Cnt > (int)TRUE_ARRAY.size())
			TRUE_ARRAY.resize(Cnt, TRUE);
	}

	vector<CBOOL> TRUE_ARRAY;
	vector<UInt8> GENO_BUFFER;
	map<int, TSelList> _Map;
} Init;




// ###########################################################
// the structure of read line
// ###########################################################

/// a class of parsing text
class CReadLine
{
public:
	/// constructor
	CReadLine()
	{
		_ReadFun = _Rho = R_NilValue;
		_ptr_line = _lines.end();
		_ifend = false; _line_no = _column_no = 0;
		_cur_char = NULL;
	}
	/// constructor
	CReadLine(SEXP vFun, SEXP vRho)
	{
		Init(vFun, vRho);
	}

	/// initialize R call
	void Init(SEXP vFun, SEXP vRho)
	{
		_ReadFun = vFun; _Rho = vRho;	
		_lines.clear(); _ptr_line = _lines.end();
		_ifend = false; _line_no = _column_no = 0;
		_cur_char = NULL;
	}

	/// read a line
	const char *ReadLine()
	{
		if (_ifend) return NULL;
		if (_ptr_line == _lines.end())
		{
			if (_PrepareBuffer())
			{
				const char *rv = *_ptr_line;
				_ptr_line ++; _line_no ++;
				return rv;
			} else
				return NULL;
		} else {
			const char *rv = *_ptr_line;
			_ptr_line ++; _line_no ++;
			return rv;
		}
	}

	/// get a string with a seperator '\t'
	void GetCell(string &buffer, bool last_column)
	{
		if (_ifend)
			throw ErrSeqArray("It is the end.");
		if (!_cur_char)
		{
			_cur_char = ReadLine();
			if (!_cur_char)
				throw ErrSeqArray("It is the end.");
			_column_no = 0;
		}

		const char *str_begin = _cur_char;
		while ((*_cur_char != '\t') && (*_cur_char != 0))
			_cur_char ++;
		const char *str_end = _cur_char;
		_column_no ++;

		// check
		if ((str_begin == str_end) && (*_cur_char == 0))
			throw ErrSeqArray("fewer columns than what expected.");
		if (last_column)
		{
			if (*_cur_char != 0)
				throw ErrSeqArray("more columns than what expected.");
			_cur_char = NULL;
		} else {
			if (*_cur_char == '\t') _cur_char ++;
		}
		buffer.assign(str_begin, str_end);
	}

	/// return true, if it is of the end
	bool IfEnd()
	{
		if (!_ifend)
		{
			if (_ptr_line == _lines.end())
				_PrepareBuffer();
		}
		return _ifend;
	}

	/// return line number
	COREARRAY_INLINE int LineNo() { return _line_no; }
	/// return column number
	COREARRAY_INLINE int ColumnNo() { return _column_no; }

protected:
	SEXP _ReadFun;  //< R call function
	SEXP _Rho;      //< R environment
	vector<const char *> _lines;               //< store returned string(s)
	vector<const char *>::iterator _ptr_line;  //< the pointer to _lines
	bool _ifend;     //< true for the end of reading
	int _line_no;    //< the index of current line
	int _column_no;  //< the index of current column
	const char *_cur_char;  //< 

	bool _PrepareBuffer()
	{
		// call ReadLine R function
		SEXP val = eval(_ReadFun, _Rho);

		// check the returned value
		int n = Rf_length(val);
		if (n > 0)
		{
			_ifend = false;
			_lines.resize(n);
			for (int i=0; i < n; i++)
				_lines[i] = CHAR(STRING_ELT(val, i));
			_ptr_line = _lines.begin();
			return true;
		} else {
			_ifend = true;
			return false;
		}
	}
};




// ###########################################################
// VCF strcture
// ###########################################################

const static int FIELD_TYPE_INT      = 1;
const static int FIELD_TYPE_FLOAT    = 2;
const static int FIELD_TYPE_FLAG     = 3;
const static int FIELD_TYPE_STRING   = 4;


/// the structure of INFO field
struct TVCF_Field_Info
{
	string name;           //< INFO ID
	int type;              //< 1: integer, 2: float, 3: flag, 4: character,
	bool import_flag;      //< true: import, false: not import
	PdSequenceX data_obj;  //< the pointer to data object
	PdSequenceX len_obj;   //< can be NULL if variable-length object
	int number;            //< according to 'Number' field, if -1: variable-length, -2: # of alleles, -3: # of genotypes
	bool used;             //< if TRUE, it has been parsed for the current line

	TVCF_Field_Info() { type = 0; data_obj = len_obj = NULL; number = 0; used = false; }

	// INFO field

	template<typename TYPE> void Check(vector<TYPE> &array, string &name, int num_allele)
	{
		CoreArray::Int32 I32;
		switch (number)
		{
			case -1:
				// variable-length
				I32 = array.size();
				CHECK(gds_AppendData(len_obj, 1, &I32, svInt32));
				break;		
			case -2:
				// # of alleles
				I32 = array.size();
				if (I32 != (num_allele-1))
				{
					throw ErrSeqArray("INFO ID '%s' should have %d value(s).",
						name.c_str(), num_allele-1);
				}
				CHECK(gds_AppendData(len_obj, 1, &I32, svInt32));
				break;		
			case -3:
				// # of genotypes
				I32 = array.size();
				if (I32 != (num_allele+1)*num_allele/2)
				{
					throw ErrSeqArray("INFO ID '%s' should have %d value(s).",
						name.c_str(), (num_allele+1)*num_allele/2);
				}
				CHECK(gds_AppendData(len_obj, 1, &I32, svInt32));
				break;		
			default:
				if (number >= 0)
				{
					if (number != (int)array.size())
					{
						throw ErrSeqArray("INFO ID '%s' should have %d value(s).",
							name.c_str(), number);
					}
				} else
					throw ErrSeqArray("Invalid value 'number' in TVCF_Field_Info.");
		}
	}

	template<typename TYPE> void Fill(vector<TYPE> &array, TYPE val)
	{
		if (number < 0)
		{
			CoreArray::Int32 I32 = 0;
			CHECK(gds_AppendData(len_obj, 1, &I32, svInt32));
		} else {
			array.clear();
			array.resize(number, val);
			CHECK(gds_AppendData(data_obj, number, &(array[0]), TdTraits<TYPE>::SVType));
		}
	}
};


/// the structure of FORMAT field
struct TVCF_Field_Format
{
	string name;           //< FORMAT ID
	int type;              //< 1: integer, 2: float, 3: flag, 4: character,
	bool import_flag;      //< true: import, false: not import
	PdSequenceX data_obj;  //< the pointer to data object
	PdSequenceX len_obj;   //< can be NULL if variable-length object
	int number;            //< according to 'Number' field, if -1: variable-length, -2: # of alleles, -3: # of genotypes
	bool used;             //< if TRUE, it has been parsed for the current line

	/// data -- Int32
	vector< vector<CoreArray::Int32> > I32ss;
	/// data -- Float32
	vector< vector<Float32> > F32ss;
	/// data -- UTF8 string
	vector< vector<string> > UTF8ss;


	TVCF_Field_Format() { type = 0; data_obj = len_obj = NULL; number = 0; used = false; }

	// FORMAT field

	template<typename TYPE>
		void Check(vector<TYPE> &array, string &name, int num_allele, const TYPE &missing)
	{
		switch (number)
		{
			case -1:
				break;

			case -2:
				// # of alleles
				if ((int)array.size() > (num_allele-1))
				{
					throw ErrSeqArray("FORMAT ID '%s' should have %d value(s).",
						name.c_str(), num_allele-1);
				} else {
					array.resize(num_allele-1, missing);
				}
				break;		

			case -3:
				// # of genotypes
				if ((int)array.size() > (num_allele+1)*num_allele/2)
				{
					throw ErrSeqArray("INFO ID '%s' should have %d value(s).",
						name.c_str(), (num_allele+1)*num_allele/2);
				} else {
					array.resize((num_allele+1)*num_allele/2, missing);
				}
				break;		

			default:
				if (number >= 0)
				{
					if ((int)array.size() > number)
					{
						throw ErrSeqArray("FORMAT ID '%s' should have %d value(s).",
							name.c_str(), number);
					} else {
						array.resize(number, missing);
					}
				} else
					throw ErrSeqArray("Invalid value 'number' in TVCF_Field_Format.");
		}
	}

	void WriteFixedLength()
	{
		if (number < 0)
			throw ErrSeqArray("Wrong call 'WriteFixedLength' in TVCF_Field_Format.");
		switch (type)
		{
			case FIELD_TYPE_INT:
				for (vector< vector<CoreArray::Int32> >::iterator it = I32ss.begin();
					it != I32ss.end(); it ++)
				{
					CHECK(gds_AppendData(data_obj, number, &((*it)[0]), svInt32));
				}
				break;

			case FIELD_TYPE_FLOAT:
				for (vector< vector<float> >::iterator it = F32ss.begin();
					it != F32ss.end(); it ++)
				{
					CHECK(gds_AppendData(data_obj, number, &((*it)[0]), svFloat32));
				}
				break;

			case FIELD_TYPE_STRING:
				for (vector< vector<string> >::iterator it = UTF8ss.begin();
					it != UTF8ss.end(); it ++)
				{
					for (int j=0; j < (int)(*it).size(); j ++)
						CHECK(gds_AppendString(data_obj, (*it)[j].c_str()));
				}
				break;

			default:
				throw ErrSeqArray("Invalid FORMAT Type.");
		}
	}

	int WriteVariableLength(int nTotalSample, vector<CoreArray::Int32> &I32s,
		vector<float> &F32s)
	{
		if (number >= 0)
			throw ErrSeqArray("Wrong call 'WriteVariableLength' in TVCF_Field_Format.");

		int nMax = 0;
		switch (type)
		{
			case FIELD_TYPE_INT:
				for (int j=0; j < nTotalSample; j++)
				{
					if (nMax < (int)I32ss[j].size())
						nMax = I32ss[j].size();
				}
				I32s.resize(nTotalSample);
				for (int i=0; i < nMax; i++)
				{
					for (int j=0; j < nTotalSample; j++)
					{
						vector<CoreArray::Int32> &B = I32ss[j];
						I32s[j] = (i < (int)B.size()) ? B[i] : NA_INTEGER;
					}
					CHECK(gds_AppendData(data_obj, nTotalSample, &(I32s[0]), svInt32));
				}
				break;

			case FIELD_TYPE_FLOAT:
				for (int j=0; j < nTotalSample; j++)
				{
					if (nMax < (int)F32ss[j].size())
						nMax = F32ss[j].size();
				}
				F32s.resize(nTotalSample);
				for (int i=0; i < nMax; i++)
				{
					for (int j=0; j < nTotalSample; j++)
					{
						vector<float> &B = F32ss[j];
						F32s[j] = (i < (int)B.size()) ? B[i] : (float)R_NaN;
					}
					CHECK(gds_AppendData(data_obj, nTotalSample, &(F32s[0]), svFloat32));
				}
				break;

			case FIELD_TYPE_STRING:
				for (int j=0; j < nTotalSample; j++)
				{
					if (nMax < (int)UTF8ss[j].size())
						nMax = UTF8ss[j].size();
				}
				for (int i=0; i < nMax; i++)
				{
					for (int j=0; j < nTotalSample; j++)
					{
						vector<string> &B = UTF8ss[j];
						CHECK(gds_AppendString(data_obj,
							(i < (int)B.size()) ? B[i].c_str() : ""));
					}
				}
				break;

			default:
				throw ErrSeqArray("Invalid FORMAT Type.");
		}
		return nMax;
	}
};



/// 
static void _MappingIndex(PdSequenceX Node, const vector<CBOOL> &sel,
	vector<int> &out_len, vector<CBOOL> &out_var_sel,
	CoreArray::Int32 &out_var_start, CoreArray::Int32 &out_var_count)
{
	if (gds_SeqDimCnt(Node) != 1)
		throw ErrSeqArray("Invalid dimension.");
	int Cnt = CHECK(gds_SeqGetCount(Node));

	if (sel.empty())
	{
		out_len.resize(Cnt);
		CoreArray::Int32 _st=0, _cnt=Cnt;
		CHECK(gds_rData(Node, &_st, &_cnt, &out_len[0], svInt32));

		out_var_start = 0;
		out_var_count = 0;
		for (vector<int>::iterator it=out_len.begin();
			it != out_len.end(); it++)
		{
			if (*it > 0) out_var_count += *it;
		}
		out_var_sel.clear();
		out_var_sel.resize(out_var_count, TRUE);

	} else {
		// check
		if ((int)sel.size() != Cnt)
			throw ErrSeqArray("Invalid dimension.");

		// find the start
		int _start = 0;
		for (; _start < (int)sel.size(); _start++)
			if (sel[_start]) break;
		// find the end
		int _end = sel.size()-1;
		for (; _end >= 0; _end --)
			if (sel[_end]) break;

		if (_end >= 0)
		{
			const int N_MAX = 16384;
			CoreArray::Int32 buffer[N_MAX];

			out_var_start = 0;
			int pos = 0;
			while (pos < _start)
			{
				int L = _start - pos;
				if (L > N_MAX) L = N_MAX;
				CHECK(gds_rData(Node, &pos, &L, buffer, svInt32));
				pos += L;
				for (int i=0; i < L; i++)
				{
					if (buffer[i] > 0)
						out_var_start += buffer[i];
				}
			}

			out_len.clear();
			out_var_sel.clear();
			while (pos <= _end)
			{
				int L = _end - pos + 1;
				if (L > N_MAX) L = N_MAX;
				CHECK(gds_rData(Node, &pos, &L, buffer, svInt32));
				for (int i=0; i < L; i++)
				{
					int LL = (buffer[i] > 0) ? buffer[i] : 0;
					if (sel[pos+i])
					{
						out_len.push_back(LL);
						for (int j=0; j < LL; j++)
							out_var_sel.push_back(TRUE);
					} else {
						for (int j=0; j < LL; j++)
							out_var_sel.push_back(FALSE);
					}
				}
				pos += L;
			}
			
			out_var_count = out_var_sel.size();
		} else {
			out_len.clear(); out_var_sel.clear();
			out_var_start = out_var_count = 0;
		}
	}
}




/// 
class TVariable_ApplyByVariant
{
public:
	enum TType {
		ctNone, ctBasic, ctGenotype, ctPhase, ctInfo, ctFormat
	};


	TVariable_ApplyByVariant()
	{
		Node = IndexNode = NULL;
		VariantSelection = NULL;
	}

	void InitObject(TType Type, const char *Path, PdGDSObj Root,
		int nVariant, CBOOL *VariantSel, int nSample, CBOOL *SampleSel)
	{
		static const char *ErrDim = "Invalid dimension of '%s'.";

		// initialize
		GDS_PATH_PREFIX_CHECK(Path);
		VarType = Type;
		Node = CHECK(gds_NodePath(Root, Path));
		SVType = gds_SeqSVType(Node);
		DimCnt = gds_SeqDimCnt(Node);

		Num_Variant = nVariant;
		TotalNum_Sample = nSample;
		VariantSelection = VariantSel;

		Num_Sample = 0;
		for (int i=0; i < nSample; i ++)
		{
			if (SampleSel[i])
				Num_Sample ++;
		}

		string Path2; // the path with '@'

		switch (Type)
		{
			case ctBasic:
				// VARIABLE: variant.id, position, allele
				if ((DimCnt != 1) || (gds_SeqGetCount(Node) != nVariant))
					throw ErrSeqArray(ErrDim, Path);
				break;

			case ctGenotype:
				// VARIABLE: genotype/data, genotype/@data
				if (DimCnt != 3)
					throw ErrSeqArray(ErrDim, Path);
				CHECK(gds_SeqGetDim(Node, DLen));
				if ((DLen[0] < nVariant) || (DLen[1] != nSample))
					throw ErrSeqArray(ErrDim, Path);

				Path2 = GDS_PATH_PREFIX(Path, '@');
				IndexNode = gds_NodePath(Root, Path2.c_str());
				if (IndexNode == NULL)
					throw ErrSeqArray("'%s' is missing!", Path2.c_str());
				if ((gds_SeqDimCnt(IndexNode) != 1) || (gds_SeqGetCount(IndexNode) != nVariant))
					throw ErrSeqArray(ErrDim, Path2.c_str());

				SelPtr[1] = SampleSel;
				Init.Check_TrueArray(DLen[2]);
				SelPtr[2] = &Init.TRUE_ARRAY[0];
				break;

			case ctPhase:
				// VARIABLE: phase/data
				if ((DimCnt != 2) && (DimCnt != 3))
					throw ErrSeqArray(ErrDim, Path);
				CHECK(gds_SeqGetDim(Node, DLen));
				if ((DLen[0] != nVariant) || (DLen[1] != nSample))
					throw ErrSeqArray(ErrDim, Path);

				SelPtr[1] = SampleSel;
				if (DimCnt > 2)
				{
					Init.Check_TrueArray(DLen[2]);
					SelPtr[2] = &Init.TRUE_ARRAY[0];  //< ToDo: check
				}
				break;

			case ctInfo:
				// VARIABLE: info/...
				if ((DimCnt!=1) && (DimCnt!=2))
					throw ErrSeqArray(ErrDim, Path);
				CHECK(gds_SeqGetDim(Node, DLen));

				Path2 = GDS_PATH_PREFIX(Path, '@');
				IndexNode = gds_NodePath(Root, Path2.c_str());
				if (IndexNode != NULL)
				{
					if ((gds_SeqDimCnt(IndexNode) != 1) || (gds_SeqGetCount(IndexNode) != nVariant))
						throw ErrSeqArray(ErrDim, Path2.c_str());
				} else {
					if (DLen[0] != nVariant)
						throw ErrSeqArray(ErrDim, Path);
				}

				if (DimCnt > 1)
				{
					Init.Check_TrueArray(DLen[1]);
					SelPtr[1] = &Init.TRUE_ARRAY[0];
				}
				break;

			case ctFormat:
				// VARIABLE: format/...
				if ((DimCnt!=2) && (DimCnt!=3))
					throw ErrSeqArray(ErrDim, Path);
				CHECK(gds_SeqGetDim(Node, DLen));

				Path2 = GDS_PATH_PREFIX(Path, '@');
				IndexNode = gds_NodePath(Root, Path2.c_str());
				if (IndexNode != NULL)
				{
					if ((gds_SeqDimCnt(IndexNode) != 1) || (gds_SeqGetCount(IndexNode) != nVariant))
						throw ErrSeqArray(ErrDim, Path2.c_str());
				} else
					throw ErrSeqArray("'%s' is missing!", Path2.c_str());

				SelPtr[1] = SampleSel;
				if (DimCnt > 2)
				{
					Init.Check_TrueArray(DLen[2]);
					SelPtr[2] = &Init.TRUE_ARRAY[0];
				}
				break;

			default:
				throw ErrSeqArray("Internal Error in 'TVariable_ApplyByVariant::InitObject'.");
		}

		_Index = 0;
		IndexCellVariant = 0;
		if (IndexNode)
		{
			CoreArray::Int32 Cnt=1;
			CHECK(gds_rData(IndexNode, &_Index, &Cnt, &NumCellVariant, svInt32));
			if (NumCellVariant < 0) NumCellVariant = 0;
		} else
			NumCellVariant = 1;
		if (!VariantSelection[0]) NextCell();

		if (Type == ctGenotype)
		{
			CellCount = Num_Sample * DLen[2];
			int SlideCnt = DLen[1] * DLen[2];
			if (SlideCnt > (int)Init.GENO_BUFFER.size())
				Init.GENO_BUFFER.resize(SlideCnt);
		}
	}

	bool NextCell()
	{
		_Index ++;
		IndexCellVariant += NumCellVariant;
		if (IndexNode)
		{
			CoreArray::Int32 Cnt=1, L;
			while ((_Index<Num_Variant) && !VariantSelection[_Index])
			{
				CHECK(gds_rData(IndexNode, &_Index, &Cnt, &L, svInt32));
				if (L > 0) IndexCellVariant += L;
				_Index ++;
			}
			if (_Index < Num_Variant)
			{
				CHECK(gds_rData(IndexNode, &_Index, &Cnt, &NumCellVariant, svInt32));
				if (NumCellVariant < 0) NumCellVariant = 0;
			} else
				NumCellVariant = 0;
		} else {
			while ((_Index<Num_Variant) && !VariantSelection[_Index])
				_Index ++;
			IndexCellVariant = _Index;
			NumCellVariant = 1;
		}
		return (_Index < Num_Variant);
	}

	void ReadGenoData(int *Base)
	{
		// the size of Init.GENO_BUFFER has been check in 'Init()'
		int SlideCnt = DLen[1]*DLen[2];

		TdIterator it;
		CHECK(gds_IterGetStart(Node, it));
		CHECK(gds_IterAdvEx(it, IndexCellVariant*SlideCnt));
		CHECK_SIZE(gds_IterRData(it, &Init.GENO_BUFFER[0], SlideCnt, svUInt8));
		UInt8 *s = &Init.GENO_BUFFER[0];
		int *p = Base;
		for (int i=0; i < DLen[1]; i++)
		{
			if (SelPtr[1][i])
			{
				for (int j=0; j < DLen[2]; j++)
					*p++ = *s++;
			} else {
				s += DLen[2];
			}
		}

		int missing = 3;

		// CellCount = Num_Sample * DLen[2] in 'NeedRData'
		for (int idx=1; idx < NumCellVariant; idx ++)
		{
			CHECK(gds_IterGetStart(Node, it));
			CHECK(gds_IterAdvEx(it, (IndexCellVariant + idx)*SlideCnt));
			CHECK_SIZE(gds_IterRData(it, &Init.GENO_BUFFER[0], SlideCnt, svUInt8));

			int shift = idx*2;
			s = &Init.GENO_BUFFER[0];
			p = Base;
			for (int i=0; i < DLen[1]; i++)
			{
				if (SelPtr[1][i])
				{
					for (int j=0; j < DLen[2]; j++)
					{
						*p |= int(*s) << shift;
						p ++; s ++;
					}
				} else {
					s += DLen[2];
				}
			}

			missing = (missing << 2) | 0x03;
		}
		for (int n=CellCount; n > 0; n--)
		{
			if (*Base == missing) *Base = NA_INTEGER;
			Base ++;
		}	
	}

	void ReadData(SEXP Val)
	{
		if (NumCellVariant <= 0) return;
		if (VarType == ctGenotype)
		{
			ReadGenoData(INTEGER(Val));
		} else {
			int st[3] = { IndexCellVariant, 0, 0 };
			DLen[0] = NumCellVariant;
			if (NumCellVariant > (int)Init.TRUE_ARRAY.size())
				Init.TRUE_ARRAY.resize(NumCellVariant);
			SelPtr[0] = &Init.TRUE_ARRAY[0];
			if (COREARRAY_SV_INTEGER(SVType))
			{
				CHECK(gds_rDataEx(Node, st, DLen, SelPtr, INTEGER(Val), svInt32));
			} else if (COREARRAY_SV_FLOAT(SVType))
			{
				CHECK(gds_rDataEx(Node, st, DLen, SelPtr, REAL(Val), svFloat64));
			} else if (COREARRAY_SV_STRING(SVType))
			{
				vector<string> buffer(CellCount);
				CHECK(gds_rDataEx(Node, st, DLen, SelPtr, &buffer[0], svStrUTF8));
				for (int i=0; i < (int)buffer.size(); i++)
					SET_STRING_ELT(Val, i, mkChar(buffer[i].c_str()));
			}
		}
	}

	SEXP NeedRData(int &nProtected)
	{
		if (NumCellVariant <= 0) return R_NilValue;

		map<int, SEXP>::iterator it = VarBuffer.find(NumCellVariant);
		if (it == VarBuffer.end())
		{
			switch (VarType)
			{
			case ctBasic:
				CellCount = 1; break;
			case ctGenotype:
				CellCount = Num_Sample * DLen[2]; break;
			case ctPhase:
				CellCount = (DimCnt>2) ? Num_Sample*DLen[2] : Num_Sample;
				break;
			case ctInfo:
				CellCount = ((DimCnt>1) ? DLen[1] : 1) * NumCellVariant;
				break;
			case ctFormat:
				CellCount = ((DimCnt>2) ? Num_Sample*DLen[2] : Num_Sample) * NumCellVariant;
				break;
			default:
				CellCount = 0;
			}

			SEXP ans = R_NilValue, dim;
			if (COREARRAY_SV_INTEGER(SVType))
			{
				char classname[32];
				classname[0] = 0;
				gds_NodeClassName(Node, classname, sizeof(classname));
				if (strcmp(classname, "dBit1") == 0)
				{
					PROTECT(ans = NEW_LOGICAL(CellCount));
				} else if (gds_Is_R_Logical(Node))
				{
					PROTECT(ans = NEW_LOGICAL(CellCount));
				} else {
					PROTECT(ans = NEW_INTEGER(CellCount));
					nProtected += gds_Set_If_R_Factor(Node, ans);
				}
				nProtected ++;
			} else if (COREARRAY_SV_FLOAT(SVType))
			{
				PROTECT(ans = NEW_NUMERIC(CellCount));
				nProtected ++;
			} else if (COREARRAY_SV_STRING(SVType))
			{
				PROTECT(ans = NEW_CHARACTER(CellCount));
				nProtected ++;
			}

			SEXP name_list, tmp;
			switch (VarType)
			{
			case ctGenotype:
				PROTECT(dim = NEW_INTEGER(2));
				INTEGER(dim)[0] = DLen[2]; INTEGER(dim)[1] = Num_Sample;
				SET_DIM(ans, dim);
				PROTECT(name_list = NEW_LIST(2));
				PROTECT(tmp = NEW_CHARACTER(2));
					SET_STRING_ELT(tmp, 0, mkChar("allele"));
					SET_STRING_ELT(tmp, 1, mkChar("sample"));
					SET_NAMES(name_list, tmp);
				SET_DIMNAMES(ans, name_list);
				nProtected += 3;
				break;

			case ctPhase:
				if (DimCnt > 2)
				{
					PROTECT(dim = NEW_INTEGER(2)); nProtected ++;
					INTEGER(dim)[0] = DLen[2]; INTEGER(dim)[1] = Num_Sample;
					SET_DIM(ans, dim);
				}
				break;

			case ctFormat:
				if (DimCnt == 2)
				{
					PROTECT(dim = NEW_INTEGER(2)); nProtected ++;
					INTEGER(dim)[0] = Num_Sample; INTEGER(dim)[1] = NumCellVariant;
					SET_DIM(ans, dim);
				} else if (DimCnt > 2)
				{
					PROTECT(dim = NEW_INTEGER(3)); nProtected ++;
					INTEGER(dim)[0] = DLen[2]; INTEGER(dim)[1] = Num_Sample;
					INTEGER(dim)[2] = NumCellVariant;
					SET_DIM(ans, dim);
				}
				break;

			default:
				break;
			}

			VarBuffer.insert(pair<int, SEXP>(NumCellVariant, ans));
			return ans;
		} else
			return it->second;
	}


	map<int, SEXP> VarBuffer;

	TType VarType;        //< 
	PdSequenceX Node;
	PdSequenceX IndexNode;
	int _Index;           //< the index of variant, starting from ZERO
	int SVType;           //< Data Type
	int DimCnt;           //< the number of dimensions
	int DLen[4];
	int Num_Variant;      //< the total number of variants
	int TotalNum_Sample;  //< the total number of samples
	int Num_Sample;       //< the number of selected samples
	CBOOL *SelPtr[3];
	CBOOL *VariantSelection;

protected:
	int IndexCellVariant;  //< 
	int NumCellVariant;    //< 
	int CellCount;         //< 
};







extern "C"
{

// ###########################################################
// the internal functions in R, "Rconnections.h"
// ###########################################################

typedef void *_Rconnection;

extern int Rconn_fgetc(_Rconnection con);
extern int Rconn_ungetc(int c, _Rconnection con);
extern _Rconnection getConnection(int n);
	



// ###########################################################
// the initial and final functions
// ###########################################################

/// initialize the package
DLLEXPORT SEXP seq_Init(SEXP lib_fn)
{
	SEXP ans = R_NilValue;

	try {
		// the file name
		const char *fn = CHAR(STRING_ELT(lib_fn, 0));
		// initialize the GDS interface
		GDSInterface::InitGDSInterface(fn);
	}
	catch (exception &E) {
		ans = mkString(E.what());
	}
	catch (const char *E) {
		ans = mkString(E);
	}
	catch (...) {
		ans = mkString("unknown error!");
	}

	return ans;
}

/// finalize the package
DLLEXPORT void seq_Done()
{
	try {
		GDSInterface::DoneGDSInterface();
	} catch (...) {};
}



// ###########################################################
// Open a GDS file
// ###########################################################

/// initialize a SeqArray file
DLLEXPORT SEXP seq_Open_Init(SEXP gdsfile)
{
	CORETRY_CALL
		TInitObject::TSelection &s = Init.Selection(gdsfile);
		s.Sample.clear();
		s.Variant.clear();
	CORECATCH_CALL
	return R_NilValue;
}

/// finalize a SeqArray file
DLLEXPORT SEXP seq_File_Done(SEXP gdsfile)
{
	CORETRY_CALL
		int gds_file_id = INTEGER(getListElement(gdsfile, "id"))[0];
		map<int, TInitObject::TSelList>::iterator it =
			Init._Map.find(gds_file_id);
		if (it != Init._Map.end())
			Init._Map.erase(it);
	CORECATCH_CALL
	return R_NilValue;
}



// ###########################################################
// Set a working space
// ###########################################################

/// push the current filter to the stack
DLLEXPORT SEXP seq_FilterPush(SEXP gdsfile)
{
	CORETRY_CALL
		int id = INTEGER(getListElement(gdsfile, "id"))[0];
		map<int, TInitObject::TSelList>::iterator it =
			Init._Map.find(id);
		if (it != Init._Map.end())
		{
			it->second.push_back(TInitObject::TSelection());
		} else
			throw ErrSeqArray("The GDS file is closed or invalid.");
	CORECATCH_CALL
	return R_NilValue;
}

/// pop up the previous filter from the stack
DLLEXPORT SEXP seq_FilterPop(SEXP gdsfile)
{
	CORETRY_CALL
		int id = INTEGER(getListElement(gdsfile, "id"))[0];
		map<int, TInitObject::TSelList>::iterator it =
			Init._Map.find(id);
		if (it != Init._Map.end())
		{
			if (it->second.size() <= 1)
				throw ErrSeqArray("No filter can be pop up.");
			it->second.pop_back();
		} else
			throw ErrSeqArray("The GDS file is closed or invalid.");
	CORECATCH_CALL
	return R_NilValue;
}

/// set a working space with selected sample id
DLLEXPORT SEXP seq_SetSpaceSample(SEXP gds, SEXP samp_sel, SEXP verbose)
{
	CORETRY_CALL

		TInitObject::TSelection &s = Init.Selection(gds);

		// the GDS root node
		PdGDSObj Root = GDS_OBJECT(getListElement(gds, "root"));
		PdSequenceX varSamp = CHECK(gds_NodePath(Root, "sample.id"));

		if (gds_SeqDimCnt(varSamp) != 1)
			throw ErrSeqArray("Invalid dimension of 'sample.id'!");
		int Count = 0;
		CHECK(gds_SeqGetDim(varSamp, &Count));

		vector<CBOOL> flag_array(Count);

		if (IS_LOGICAL(samp_sel))
		{
			// a logical vector for selected samples
			if (Rf_length(samp_sel) != Count)
				throw ErrSeqArray("Invalid length of 'samp.sel'.");
			// set selection
			int *base = LOGICAL(samp_sel);
			for (int i=0; i < Count; i++)
				flag_array[i] = (base[i] == TRUE);

		} else if (IS_INTEGER(samp_sel))
		{
			// initialize
			set<int> set_id;
			set_id.insert(INTEGER(samp_sel), &INTEGER(samp_sel)[Rf_length(samp_sel)]);
			// sample id
			vector<int> sample_id(Count);
			CoreArray::Int32 _st=0, _cnt=Count;
			CHECK(gds_rData(varSamp, &_st, &_cnt, &sample_id[0], svInt32));
			// set selection
			for (int i=0; i < Count; i++)
			{
				set<int>::iterator it = set_id.find(sample_id[i]);
				flag_array[i] = (it != set_id.end());
			}

		} else if (IS_NUMERIC(samp_sel))
		{
			// initialize
			set<double> set_id;
			set_id.insert(REAL(samp_sel), &REAL(samp_sel)[Rf_length(samp_sel)]);
			// sample id
			vector<double> sample_id(Count);
			CoreArray::Int32 _st=0, _cnt=Count;
			CHECK(gds_rData(varSamp, &_st, &_cnt, &sample_id[0], svFloat64));
			// set selection
			for (int i=0; i < Count; i++)
			{
				set<double>::iterator it = set_id.find(sample_id[i]);
				flag_array[i] = (it != set_id.end());
			}

		} else if (IS_CHARACTER(samp_sel))
		{
			// initialize
			set<string> set_id;
			for (int i=0; i < Rf_length(samp_sel); i++)
				set_id.insert(string(CHAR(STRING_ELT(samp_sel, i))));
			// sample id
			vector<string> sample_id(Count);
			CoreArray::Int32 _st=0, _cnt=Count;
			CHECK(gds_rData(varSamp, &_st, &_cnt, &sample_id[0], svStrUTF8));
			// set selection
			for (int i=0; i < Count; i++)
			{
				set<string>::iterator it = set_id.find(sample_id[i]);
				flag_array[i] = (it != set_id.end());
			}
		
		} else if (isNull(samp_sel))
		{
			flag_array.clear();
		} else
			throw ErrSeqArray("Invalid type of 'samp_sel'.");

		int n = 0;
		for (vector<CBOOL>::iterator it=flag_array.begin();
			it != flag_array.end(); it ++)
		{
			if (*it != 0) n ++;
		}
		if (isNull(samp_sel)) n = Count;
		if (n > 0)
		{
			s.Sample = flag_array;
			if (LOGICAL(verbose)[0] == TRUE)
				Rprintf("# of selected samples: %d\n", n);
		} else
			throw ErrSeqArray("No sample selected!");

	CORECATCH_CALL

	// output
	return(R_NilValue);
}


/// set a working space with selected variant id
DLLEXPORT SEXP seq_SetSpaceVariant(SEXP gds, SEXP var_sel, SEXP verbose)
{
	CORETRY_CALL

		TInitObject::TSelection &s = Init.Selection(gds);

		// the GDS root node
		PdGDSObj Root = GDS_OBJECT(getListElement(gds, "root"));
		PdSequenceX varVariant = CHECK(gds_NodePath(Root, "variant.id"));

		if (gds_SeqDimCnt(varVariant) != 1)
			throw ErrSeqArray("Invalid dimension of 'variant.id'!");
		int Count = 0;
		CHECK(gds_SeqGetDim(varVariant, &Count));

		vector<CBOOL> flag_array(Count);

		if (IS_LOGICAL(var_sel))
		{
			// a logical vector for selected samples
			if (Rf_length(var_sel) != Count)
				throw ErrSeqArray("Invalid length of 'variant.sel'.");
			// set selection
			int *base = LOGICAL(var_sel);
			for (int i=0; i < Count; i++)
				flag_array[i] = (base[i] == TRUE);

		} else if (IS_INTEGER(var_sel))
		{
			// initialize
			set<int> set_id;
			set_id.insert(INTEGER(var_sel), &INTEGER(var_sel)[Rf_length(var_sel)]);
			// sample id
			vector<int> var_id(Count);
			CoreArray::Int32 _st=0, _cnt=Count;
			CHECK(gds_rData(varVariant, &_st, &_cnt, &var_id[0], svInt32));
			// set selection
			for (int i=0; i < Count; i++)
			{
				set<int>::iterator it = set_id.find(var_id[i]);
				flag_array[i] = (it != set_id.end());
			}

		} else if (IS_NUMERIC(var_sel))
		{
			// initialize
			set<double> set_id;
			set_id.insert(REAL(var_sel), &REAL(var_sel)[Rf_length(var_sel)]);
			// variant id
			vector<double> variant_id(Count);
			CoreArray::Int32 _st=0, _cnt=Count;
			CHECK(gds_rData(varVariant, &_st, &_cnt, &variant_id[0], svFloat64));
			// set selection
			for (int i=0; i < Count; i++)
			{
				set<double>::iterator it = set_id.find(variant_id[i]);
				flag_array[i] = (it != set_id.end());
			}

		} else if (IS_CHARACTER(var_sel))
		{
			// initialize
			set<string> set_id;
			for (int i=0; i < Rf_length(var_sel); i++)
				set_id.insert(string(CHAR(STRING_ELT(var_sel, i))));
			// sample id
			vector<string> variant_id(Count);
			CoreArray::Int32 _st=0, _cnt=Count;
			CHECK(gds_rData(varVariant, &_st, &_cnt, &variant_id[0], svStrUTF8));
			// set selection
			for (int i=0; i < Count; i++)
			{
				set<string>::iterator it = set_id.find(variant_id[i]);
				flag_array[i] = (it != set_id.end());
			}

		} else if (isNull(var_sel))
		{
			flag_array.clear();
		} else
			throw ErrSeqArray("Invalid type of 'samp_sel'.");

		int n = 0;
		for (vector<CBOOL>::iterator it=flag_array.begin();
			it != flag_array.end(); it ++)
		{
			if (*it != 0) n ++;
		}
		if (isNull(var_sel)) n = Count;
		if (n > 0)
		{
			s.Variant = flag_array;
			if (LOGICAL(verbose)[0] == TRUE)
				Rprintf("# of selected variants: %d\n", n);
		} else
			throw ErrSeqArray("No variant selected!");

	CORECATCH_CALL

	// output
	return(R_NilValue);
}


static void CLEAR_SEL_VALUE(int num, vector<CBOOL>::iterator &it)
{
	while (num > 0)
	{
		if (*it != 0) { num --; *it = FALSE; }
		it ++;
	}
}
static void SKIP_SEL(int num, vector<CBOOL>::iterator &it)
{
	while (num > 0)
	{
		if (*it != 0) num --;
		it ++;
	}
}

/// split the selected variants according to multiple processes
DLLEXPORT SEXP seq_SplitSelectedVariant(SEXP gdsfile, SEXP Index, SEXP n_process)
{
	// selection object
	TInitObject::TSelection &s = Init.Selection(gdsfile);

	// the index process starting from 1
	int Process_Index = INTEGER(AS_INTEGER(Index))[0] - 1;
	int Num_Process = INTEGER(AS_INTEGER(n_process))[0];

	// the total number of selected variants
	vector<CBOOL>::iterator it;
	int N_Variant = 0;
	for (it=s.Variant.begin(); it != s.Variant.end();)
	{
		if (*it != 0) N_Variant ++;
		it ++;
	}
	if (N_Variant <= 0) error("No variant!");

	// split a list
	vector<int> split(Num_Process);
	double avg = (double)N_Variant / Num_Process;
	double start = 0;
	for (int i=0; i < Num_Process; i++)
	{
		start += avg;
		split[i] = (int)(start + 0.5);
	}

	// ***************************************************
	it = s.Variant.begin();
	int st = 0;
	for (int i=0; i < Process_Index; i++)
	{
		CLEAR_SEL_VALUE(split[i] - st, it);
		st = split[i];
	}
	int ans_n = split[Process_Index] - st;
	SKIP_SEL(ans_n, it);
	st = split[Process_Index];
	for (int i=Process_Index+1; i < Num_Process; i++)
	{
		CLEAR_SEL_VALUE(split[i] - st, it);
		st = split[i];
	}

	// ***************************************************
	// output
	SEXP rv = NEW_INTEGER(1);
	PROTECT(rv);
	INTEGER(rv)[0] = ans_n;
	UNPROTECT(1);

	return(rv);
}


/// split the selected samples according to multiple processes
DLLEXPORT SEXP seq_SplitSelectedSample(SEXP gdsfile, SEXP Index, SEXP n_process)
{
	// selection object
	TInitObject::TSelection &s = Init.Selection(gdsfile);

	// the index process starting from 1
	int Process_Index = INTEGER(AS_INTEGER(Index))[0] - 1;
	int Num_Process = INTEGER(AS_INTEGER(n_process))[0];

	// the total number of selected samples
	vector<CBOOL>::iterator it;
	int N_Sample = 0;
	for (it=s.Sample.begin(); it != s.Sample.end();)
	{
		if (*it != 0) N_Sample ++;
		it ++;
	}
	if (N_Sample <= 0) error("No sample!");

	// split a list
	vector<int> split(Num_Process);
	double avg = (double)N_Sample / Num_Process;
	double start = 0;
	for (int i=0; i < Num_Process; i++)
	{
		start += avg;
		split[i] = (int)(start + 0.5);
	}

	// ***************************************************
	it = s.Sample.begin();
	int st = 0;
	for (int i=0; i < Process_Index; i++)
	{
		CLEAR_SEL_VALUE(split[i] - st, it);
		st = split[i];
	}
	int ans_n = split[Process_Index] - st;
	SKIP_SEL(ans_n, it);
	st = split[Process_Index];
	for (int i=Process_Index+1; i < Num_Process; i++)
	{
		CLEAR_SEL_VALUE(split[i] - st, it);
		st = split[i];
	}

	// ***************************************************
	// output
	SEXP rv = NEW_INTEGER(1);
	PROTECT(rv);
	INTEGER(rv)[0] = ans_n;
	UNPROTECT(1);

	return(rv);
}


/// set a working space with selected variant id
DLLEXPORT SEXP seq_GetSpace(SEXP gdsfile)
{
	SEXP rv_ans = R_NilValue;
	CORETRY_CALL

		TInitObject::TSelection &s = Init.Selection(gdsfile);

		// the GDS root node
		PdGDSObj Root = GDS_OBJECT(getListElement(gdsfile, "root"));
		PdSequenceX varSample = CHECK(gds_NodePath(Root, "sample.id"));
		PdSequenceX varVariant = CHECK(gds_NodePath(Root, "variant.id"));

		int nProtected = 0;
		SEXP tmp;

		PROTECT(rv_ans = NEW_LIST(2));
		nProtected ++;

		if (s.Sample.empty())
		{
			int L = CHECK(gds_SeqGetCount(varSample));
			PROTECT(tmp = NEW_LOGICAL(L));
			nProtected ++;
			for (int i=0; i < L; i++) LOGICAL(tmp)[i] = TRUE;
		} else {
			PROTECT(tmp = NEW_LOGICAL(s.Sample.size()));
			nProtected ++;
			for (int i=0; i < (int)s.Sample.size(); i++)
				LOGICAL(tmp)[i] = (s.Sample[i] != 0);
		}
		SET_ELEMENT(rv_ans, 0, tmp);

		if (s.Variant.empty())
		{
			int L = CHECK(gds_SeqGetCount(varVariant));
			PROTECT(tmp = NEW_LOGICAL(L));
			nProtected ++;
			for (int i=0; i < L; i++) LOGICAL(tmp)[i] = TRUE;
		} else {
			PROTECT(tmp = NEW_LOGICAL(s.Variant.size()));
			nProtected ++;
			for (int i=0; i < (int)s.Variant.size(); i++)
				LOGICAL(tmp)[i] = (s.Variant[i] != 0);
		}
		SET_ELEMENT(rv_ans, 1, tmp);

		PROTECT(tmp = NEW_CHARACTER(2));
		nProtected ++;
			SET_STRING_ELT(tmp, 0, mkChar("sample.sel"));
			SET_STRING_ELT(tmp, 1, mkChar("variant.sel"));
		SET_NAMES(rv_ans, tmp);

		UNPROTECT(nProtected);

	CORECATCH_CALL

	// output
	return(rv_ans);
}


/// set a working space with selected variant id
DLLEXPORT SEXP seq_VarSummary(SEXP gdsfile, SEXP varname)
{
	SEXP rv_ans = R_NilValue;

	CORETRY_CALL

		// the selection
		TInitObject::TSelection &Sel = Init.Selection(gdsfile);
		// the GDS root node
		PdGDSObj Root = GDS_OBJECT(getListElement(gdsfile, "root"));
		// the variable name
		string vn = CHAR(STRING_ELT(varname, 0));

		if ((vn=="genotype") || (vn=="phase"))
		{
			PdGDSObj vSample  = CHECK(gds_NodePath(Root, "sample.id"));
			PdGDSObj vVariant = CHECK(gds_NodePath(Root, "variant.id"));
			PdGDSObj vGeno = gds_NodePath(Root, "genotype/data");
			if (vGeno == NULL)
			{
				vGeno = gds_NodePath(Root, "genotype/~data");
				if (vGeno == NULL)
					throw ErrSeqArray("There is no 'genotype/data' or 'genotype/~data'.");
			}

			PROTECT(rv_ans = NEW_LIST(2));
				SEXP I32, S32;

				PROTECT(I32 = NEW_INTEGER(3));
				SET_ELEMENT(rv_ans, 0, I32);
				int Buf[256];
				CHECK(gds_SeqGetDim(vGeno, Buf));
				INTEGER(I32)[0] = Buf[2];
				INTEGER(I32)[1] = gds_SeqGetCount(vSample);
				INTEGER(I32)[2] = gds_SeqGetCount(vVariant);

				PROTECT(S32 = NEW_INTEGER(2));
				SET_ELEMENT(rv_ans, 1, S32);
				if (!Sel.Sample.empty())
				{
					int &n = INTEGER(S32)[0]; n = 0;
					vector<CBOOL>::iterator it;
					for (it=Sel.Sample.begin(); it != Sel.Sample.end(); it ++)
						if (*it) n ++;
				} else
					INTEGER(S32)[0] = INTEGER(I32)[1];
				if (!Sel.Variant.empty())
				{
					int &n = INTEGER(S32)[1]; n = 0;
					vector<CBOOL>::iterator it;
					for (it=Sel.Variant.begin(); it != Sel.Variant.end(); it ++)
						if (*it) n ++;
				} else
					INTEGER(S32)[1] = INTEGER(I32)[2];

			SEXP tmp;
			PROTECT(tmp = NEW_CHARACTER(2));
				SET_STRING_ELT(tmp, 0, mkChar("dim"));
				SET_STRING_ELT(tmp, 1, mkChar("seldim"));
				SET_NAMES(rv_ans, tmp);
			UNPROTECT(4);
		}

	CORECATCH_CALL

	// output
	return(rv_ans);
}


/// delete the variables
DLLEXPORT SEXP seq_Delete(SEXP gds, SEXP info, SEXP format)
{
	CORETRY_CALL

		// the GDS root node
		PdGDSObj Root = GDS_OBJECT(getListElement(gds, "root"));

		// check 
		for (int i=0; i < Rf_length(info); i++)
		{
			const char *n = CHAR(STRING_ELT(info, i));
			GDS_VARIABLE_NAME_CHECK(n);
		}
		for (int i=0; i < Rf_length(format); i++)
		{
			const char *n = CHAR(STRING_ELT(format, i));
			GDS_VARIABLE_NAME_CHECK(n);
		}

		// delete
		for (int i=0; i < Rf_length(info); i++)
		{
			const char *nm = CHAR(STRING_ELT(info, i));
			PdGDSObj N;
			N = CHECK(gds_NodePath(Root, (string("annotation/info/") + nm).c_str()));
			CHECK(gds_NodeDelete(N));
			N = gds_NodePath(Root, (string("annotation/info/@") + nm).c_str());
			if (N != NULL)
				CHECK(gds_NodeDelete(N));
		}
		for (int i=0; i < Rf_length(format); i++)
		{
			const char *nm = CHAR(STRING_ELT(format, i));
			PdGDSObj N;

			N = CHECK(gds_NodePath(Root, (string("annotation/format/") + nm + "/data").c_str()));
			CHECK(gds_NodeDelete(N));
			N = CHECK(gds_NodePath(Root, (string("annotation/format/") + nm + "/@data").c_str()));
			CHECK(gds_NodeDelete(N));
			N = gds_NodePath(Root, (string("annotation/format/") + nm + "/~data").c_str());
			if (N != NULL)
				CHECK(gds_NodeDelete(N));
			N = CHECK(gds_NodePath(Root, (string("annotation/format/") + nm).c_str()));
			CHECK(gds_NodeDelete(N));
		}

	CORECATCH_CALL

	// output
	return(R_NilValue);
}


// ###########################################################
// Get data from a working space
// ###########################################################

static SEXP VarLogical(PdGDSObj Node, SEXP Array)
{
	char classname[32];
	classname[0] = 0;
	gds_NodeClassName(Node, classname, sizeof(classname));
	if (strcmp(classname, "dBit1") == 0)
	{
		PROTECT(Array);
		Array = AS_LOGICAL(Array);
		UNPROTECT(1);
	}
	return Array;
}

/// Get data from a working space
DLLEXPORT SEXP seq_GetData(SEXP gdsfile, SEXP var_name)
{
	SEXP rv_ans = R_NilValue, tmp;
	CORETRY_CALL

		// the selection
		TInitObject::TSelection &Sel = Init.Selection(gdsfile);
		// the GDS root node
		PdGDSObj Root = GDS_OBJECT(getListElement(gdsfile, "root"));

		// 
		CBOOL *SelPtr[256];
		int DimCnt, DStart[256], DLen[256];

		// the path of GDS variable
		const char *s = CHAR(STRING_ELT(var_name, 0));
		if (strcmp(s, "sample.id") == 0)
		{
			// ***********************************************************
			// sample.id

			PdSequenceX N = CHECK(gds_NodePath(Root, s));
			DimCnt = gds_SeqDimCnt(N);
			if (DimCnt != 1)
				throw ErrSeqArray("Invalid dimension of 'sample.id'.");
			if (Sel.Sample.empty())
			{
				rv_ans = _(gds_Read_SEXP(N, NULL, NULL, NULL));
			} else {
				CHECK(gds_SeqGetDim(N, DLen));
				if ((int)Sel.Sample.size() != DLen[0])
					throw ErrSeqArray("Invalid dimension of 'sample.id'.");
				SelPtr[0] = &Sel.Sample[0];
				rv_ans = _(gds_Read_SEXP(N, NULL, NULL, &SelPtr[0]));
			}

		} else if ( (strcmp(s, "variant.id")==0) || (strcmp(s, "position")==0) ||
			(strcmp(s, "chromosome")==0) || (strcmp(s, "allele")==0) ||
			(strcmp(s, "annotation/id")==0) || (strcmp(s, "annotation/qual")==0) ||
			(strcmp(s, "annotation/filter")==0) )
		{
			// ***********************************************************
			// variant.id, position, chromosome, allele, annotation/id
			// annotation/qual, annotation/filter

			PdSequenceX N = CHECK(gds_NodePath(Root, s));
			DimCnt = gds_SeqDimCnt(N);
			if (DimCnt != 1)
				throw ErrSeqArray("Invalid dimension of '%s'.", s);
			if (Sel.Variant.empty())
			{
				rv_ans = _(gds_Read_SEXP(N, NULL, NULL, NULL));
			} else {
				CHECK(gds_SeqGetDim(N, DLen));
				if ((int)Sel.Variant.size() != DLen[0])
					throw ErrSeqArray("Invalid dimension of '%s'.", s);
				SelPtr[0] = &Sel.Variant[0];
				rv_ans = _(gds_Read_SEXP(N, NULL, NULL, &SelPtr[0]));
			}

		} else if (strcmp(s, "phase") == 0)
		{
			// *******************************************************
			// phase/

			PdSequenceX N = CHECK(gds_NodePath(Root, "phase/data"));
			DimCnt = gds_SeqDimCnt(N);
			if ((DimCnt != 2) && (DimCnt != 3))
				throw ErrSeqArray("Invalid dimension of '%s'.", s);
			if (!Sel.Sample.empty() || !Sel.Variant.empty())
			{
				CHECK(gds_SeqGetDim(N, DLen));

				if (Sel.Variant.empty())
					Sel.Variant.resize(DLen[0], TRUE);
				else if ((int)Sel.Variant.size() != DLen[0])
					throw ErrSeqArray("Invalid dimension of '%s'.", s);

				if (Sel.Sample.empty())
					Sel.Sample.resize(DLen[1], TRUE);
				else if ((int)Sel.Sample.size() != DLen[1])
					throw ErrSeqArray("Invalid dimension of '%s'.", s);

				SelPtr[0] = &Sel.Variant[0];
				SelPtr[1] = &Sel.Sample[0];
				if (DimCnt == 3)
				{
					Init.Check_TrueArray(DLen[2]);
					SelPtr[2] = &Init.TRUE_ARRAY[0];
				}

				rv_ans = _(gds_Read_SEXP(N, NULL, NULL, &SelPtr[0]));
			} else {
				rv_ans = _(gds_Read_SEXP(N, NULL, NULL, NULL));
			}

		} else if (strcmp(s, "genotype") == 0)
		{
			// *******************************************************
			// genotypic data

			// init selection
			if (Sel.Sample.empty())
			{
				PdSequenceX N = CHECK(gds_NodePath(Root, "sample.id"));
				int Cnt = gds_SeqGetCount(N);
				if (Cnt < 0) throw ErrSeqArray("Invalid dimension of 'sample.id'.");
				Sel.Sample.resize(Cnt, TRUE);
			}
			if (Sel.Variant.empty())
			{
				PdSequenceX N = CHECK(gds_NodePath(Root, "variant.id"));
				int Cnt = gds_SeqGetCount(N);
				if (Cnt < 0) throw ErrSeqArray("Invalid dimension of 'variant.id'.");
				Sel.Variant.resize(Cnt, TRUE);
			}

			// the number of selected variants
			int nVariant = 0;
			for (vector<CBOOL>::iterator it = Sel.Variant.begin();
				it != Sel.Variant.end(); it ++)
			{
				if (*it) nVariant ++;
			}
			if (nVariant > 0)
			{
				// initialize the GDS Node list
				TVariable_ApplyByVariant NodeVar;
				NodeVar.InitObject(TVariable_ApplyByVariant::ctGenotype,
					"genotype/data", Root, Sel.Variant.size(),
					&Sel.Variant[0], Sel.Sample.size(), &Sel.Sample[0]);

				// the number of calling PROTECT
				int SIZE = NodeVar.Num_Sample * NodeVar.DLen[2];
				PROTECT(rv_ans = NEW_INTEGER(nVariant * SIZE));
				PROTECT(tmp = NEW_INTEGER(3));
					INTEGER(tmp)[0] = NodeVar.DLen[2];
					INTEGER(tmp)[1] = NodeVar.Num_Sample;
					INTEGER(tmp)[2] = nVariant;
				SET_DIM(rv_ans, tmp);
				SEXP name_list;
				PROTECT(name_list = NEW_LIST(3));
				PROTECT(tmp = NEW_CHARACTER(3));
					SET_STRING_ELT(tmp, 0, mkChar("allele"));
					SET_STRING_ELT(tmp, 1, mkChar("sample"));
					SET_STRING_ELT(tmp, 2, mkChar("variant"));
					SET_NAMES(name_list, tmp);
				SET_DIMNAMES(rv_ans, name_list);

				int *base = INTEGER(rv_ans);
				do {
					NodeVar.ReadGenoData(base);
					base += SIZE;
				} while (NodeVar.NextCell());

				// finally
				UNPROTECT(4);
			}

		} else if (strncmp(s, "annotation/info/", 16) == 0)
		{
			GDS_PATH_PREFIX_CHECK(s);
			PdSequenceX N = CHECK(gds_NodePath(Root, s));
			DimCnt = gds_SeqDimCnt(N);
			if ((DimCnt!=1) && (DimCnt!=2))
				throw ErrSeqArray("Invalid dimension of '%s'.", s);

			string path_ex = GDS_PATH_PREFIX(s, '@');
			PdSequenceX N_idx = gds_NodePath(Root, path_ex.c_str());
			if (N_idx == NULL)
			{
				// no index
				if (!Sel.Variant.empty())
				{
					CHECK(gds_SeqGetDim(N, DLen));
					SelPtr[0] = &Sel.Variant[0];
					if (DimCnt == 2)
					{
						Init.Check_TrueArray(DLen[1]);
						SelPtr[1] = &Init.TRUE_ARRAY[0];
					}
	
					rv_ans = _(gds_Read_SEXP(N, NULL, NULL, &SelPtr[0]));
				} else
					rv_ans = _(gds_Read_SEXP(N, NULL, NULL, NULL));

				rv_ans = VarLogical(N, rv_ans);

			} else {
				// with index
				if (!Sel.Variant.empty())
				{
					memset(DStart, 0, sizeof(DStart));
					CHECK(gds_SeqGetDim(N, DLen));

					vector<int> len;
					vector<CBOOL> var_sel;
					_MappingIndex(N_idx, Sel.Variant, len, var_sel, DStart[0], DLen[0]);

					SelPtr[0] = &var_sel[0];
					if (DimCnt == 2)
					{
						Init.Check_TrueArray(DLen[1]);
						SelPtr[1] = &Init.TRUE_ARRAY[0];
					}

					PROTECT(rv_ans = NEW_LIST(2));
						SEXP I32;
						PROTECT(I32 = NEW_INTEGER(len.size()));
						int *base = INTEGER(I32);
						for (int i=0; i < (int)len.size(); i++)
							base[i] = len[i];
						SET_ELEMENT(rv_ans, 0, I32);
						SET_ELEMENT(rv_ans, 1,
							VarLogical(N, _(gds_Read_SEXP(N, DStart, DLen, &SelPtr[0]))));
					PROTECT(tmp = NEW_CHARACTER(2));
						SET_STRING_ELT(tmp, 0, mkChar("length"));
						SET_STRING_ELT(tmp, 1, mkChar("data"));
						SET_NAMES(rv_ans, tmp);
					UNPROTECT(3);

				} else {
					PROTECT(rv_ans = NEW_LIST(2));
						SET_ELEMENT(rv_ans, 0,
							_(gds_Read_SEXP(N_idx, NULL, NULL, NULL)));
						SET_ELEMENT(rv_ans, 1,
							VarLogical(N, _(gds_Read_SEXP(N, NULL, NULL, NULL))));
					PROTECT(tmp = NEW_CHARACTER(2));
						SET_STRING_ELT(tmp, 0, mkChar("length"));
						SET_STRING_ELT(tmp, 1, mkChar("data"));
						SET_NAMES(rv_ans, tmp);
					UNPROTECT(2);
				}
			}

		} else if (strncmp(s, "annotation/format/", 18) == 0)
		{
			GDS_PATH_PREFIX_CHECK(s);
			PdSequenceX N =
				CHECK(gds_NodePath(Root, string(string(s)+"/data").c_str()));
			PdSequenceX N_idx =
				CHECK(gds_NodePath(Root, string(string(s)+"/@data").c_str()));

			DimCnt = gds_SeqDimCnt(N);
			if ((DimCnt!=2) && (DimCnt!=3))
				throw ErrSeqArray("Invalid dimension of '%s'.", s);
			memset(DStart, 0, sizeof(DStart));
			CHECK(gds_SeqGetDim(N, DLen));

			if (Sel.Sample.empty())
				Sel.Sample.resize(DLen[1], TRUE);
			if (Sel.Variant.empty())
				Sel.Variant.resize(gds_SeqGetCount(N_idx), TRUE);

			vector<int> len;
			vector<CBOOL> var_sel;
			_MappingIndex(N_idx, Sel.Variant, len, var_sel, DStart[0], DLen[0]);

			SelPtr[0] = &var_sel[0];
			SelPtr[1] = &Sel.Sample[0];
			if (DimCnt == 3)
			{
				Init.Check_TrueArray(DLen[2]);
				SelPtr[2] = &Init.TRUE_ARRAY[0];
			}

			PROTECT(rv_ans = NEW_LIST(2));
				SEXP I32;
				PROTECT(I32 = NEW_INTEGER(len.size()));
				int *base = INTEGER(I32);
				for (int i=0; i < (int)len.size(); i++)
					base[i] = len[i];
				SET_ELEMENT(rv_ans, 0, I32);
				SEXP DAT = _(gds_Read_SEXP(N, DStart, DLen, &SelPtr[0]));
				SET_ELEMENT(rv_ans, 1, DAT);
			PROTECT(tmp = NEW_CHARACTER(2));
				SET_STRING_ELT(tmp, 0, mkChar("length"));
				SET_STRING_ELT(tmp, 1, mkChar("data"));
				SET_NAMES(rv_ans, tmp);

				if (Rf_length(DAT) > 0)
				{
					SEXP name_list;
					PROTECT(name_list = NEW_LIST(DimCnt));
					PROTECT(tmp = NEW_CHARACTER(DimCnt));
						SET_STRING_ELT(tmp, 0, mkChar("sample"));
						SET_STRING_ELT(tmp, 1, mkChar("variant"));
						SET_NAMES(name_list, tmp);
					SET_DIMNAMES(VECTOR_ELT(rv_ans, 1), name_list);
					UNPROTECT(5);
				} else {
					UNPROTECT(3);
				}

		} else {
			throw ErrSeqArray("'%s' is not a standard variable name, and the standard format:\n"
				"\tsample.id, variant.id, position, chromosome, allele, "
				"annotation/id, annotation/qual, annotation/filter\n"
				"\tannotation/info/VARIABLE_NAME, annotation/format/VARIABLE_NAME", s);
		}

	CORECATCH_CALL

	// output
	return(rv_ans);
}



// ###########################################################
// Apply functions over margins on a working space
// ###########################################################

/// Apply functions over margins on a working space
DLLEXPORT SEXP seq_Apply_Variant(SEXP gdsfile, SEXP var_name,
	SEXP FUN, SEXP as_is, SEXP var_index, SEXP rho)
{
	SEXP rv_ans = R_NilValue;
	CORETRY_CALL

		// the selection
		TInitObject::TSelection &Sel = Init.Selection(gdsfile);
		// the GDS root node
		PdGDSObj Root = GDS_OBJECT(getListElement(gdsfile, "root"));

		// init selection
		if (Sel.Sample.empty())
		{
			PdSequenceX N = CHECK(gds_NodePath(Root, "sample.id"));
			int Cnt = gds_SeqGetCount(N);
			if (Cnt < 0) throw ErrSeqArray("Invalid dimension of 'sample.id'.");
			Sel.Sample.resize(Cnt, TRUE);
		}
		if (Sel.Variant.empty())
		{
			PdSequenceX N = CHECK(gds_NodePath(Root, "variant.id"));
			int Cnt = gds_SeqGetCount(N);
			if (Cnt < 0) throw ErrSeqArray("Invalid dimension of 'variant.id'.");
			Sel.Variant.resize(Cnt, TRUE);
		}

		// the number of calling PROTECT
		int nProtected = 0;
		// the number of selected variants
		int nVariant = 0;
		for (vector<CBOOL>::iterator it = Sel.Variant.begin();
			it != Sel.Variant.end(); it ++)
		{
			if (*it) nVariant ++;
		}
		if (nVariant <= 0)
			throw ErrSeqArray("There is no selected variant.");


		// ***************************************************************
		// initialize the GDS Node list

		vector<TVariable_ApplyByVariant> NodeList(Rf_length(var_name));
		vector<TVariable_ApplyByVariant>::iterator it;

		// for - loop
		for (int i=0; i < Rf_length(var_name); i++)
		{
			// the path of GDS variable
			string s = CHAR(STRING_ELT(var_name, i));
			TVariable_ApplyByVariant::TType VarType;

			if ( s=="variant.id" || s=="position" || s=="chromosome" ||
				s=="allele" || s=="annotation/id" || s=="annotation/qual" ||
				s=="annotation/filter" )
			{
				// ***********************************************************
				// variant.id, position, chromosome, allele, annotation/id
				// annotation/qual, annotation/filter
				VarType = TVariable_ApplyByVariant::ctBasic;
			} else if (s == "genotype")
			{
				VarType = TVariable_ApplyByVariant::ctGenotype;
				s.append("/data");
			} else if (s == "phase")
			{
				// *******************************************************
				// phase/
				VarType = TVariable_ApplyByVariant::ctPhase;
				s.append("/data");
			} else if (strncmp(s.c_str(), "annotation/info/", 16) == 0)
			{
				VarType = TVariable_ApplyByVariant::ctInfo;
			} else if (strncmp(s.c_str(), "annotation/format/", 18) == 0)
			{
				VarType = TVariable_ApplyByVariant::ctFormat;
				s.append("/data");
			} else {
				throw ErrSeqArray("'%s' is not a standard variable name, and the standard format:\n"
					"\tvariant.id, position, chromosome, allele, "
					"annotation/id, annotation/qual, annotation/filter\n"
					"\tannotation/info/VARIABLE_NAME', annotation/format/VARIABLE_NAME",
					s.c_str());
			}

			NodeList[i].InitObject(VarType, s.c_str(), Root, Sel.Variant.size(),
				&Sel.Variant[0], Sel.Sample.size(), &Sel.Sample[0]);
		}

		// ***********************************************************
		// as.is
		//     0: integer, 1: double, 2: character, 3: list, other: NULL
		int DatType;
		const char *as = CHAR(STRING_ELT(as_is, 0));
		if (strcmp(as, "integer") == 0)
			DatType = 0;
		else if (strcmp(as, "double") == 0)
			DatType = 1;
		else if (strcmp(as, "character") == 0)
			DatType = 2;
		else if (strcmp(as, "list") == 0)
			DatType = 3;
		else if (strcmp(as, "none") == 0)
			DatType = -1;
		else
			throw ErrSeqArray("'as.is' is not valid!");

		// init return values
		// int DatType;  //< 0: integer, 1: double, 2: character, 3: list, other: NULL
		switch (DatType)
		{
		case 0:
			PROTECT(rv_ans = NEW_INTEGER(nVariant)); nProtected ++;
			break;
		case 1:
			PROTECT(rv_ans = NEW_NUMERIC(nVariant)); nProtected ++;
			break;
		case 2:
			PROTECT(rv_ans = NEW_CHARACTER(nVariant)); nProtected ++;
			break;
		case 3:
			PROTECT(rv_ans = NEW_LIST(nVariant)); nProtected ++;
			break;
		default:
			rv_ans = R_NilValue;
		}

		// ***********************************************************
		// rho
		if (!isEnvironment(rho))
			throw ErrSeqArray("'rho' should be an environment");


		// ***************************************************************
		// initialize calling

		SEXP R_call_param = R_NilValue;
		if (NodeList.size() > 1)
		{
			PROTECT(R_call_param = NEW_LIST(NodeList.size()));
			nProtected ++;
			// set name to R_call_param
			SET_NAMES(R_call_param, GET_NAMES(var_name));
		}

		// 1 -- none, 2 -- relative, 3 -- absolute
		int VarIdx = INTEGER(var_index)[0];

		SEXP R_fcall;
		SEXP R_Index = NULL;
		if (VarIdx > 1)
		{
			PROTECT(R_Index = NEW_INTEGER(1));
			nProtected ++;
			PROTECT(R_fcall = LCONS(FUN, LCONS(R_Index,
				LCONS(R_call_param, LCONS(R_DotsSymbol, R_NilValue)))));
			nProtected ++;
		} else {
			PROTECT(R_fcall = LCONS(FUN,
				LCONS(R_call_param, LCONS(R_DotsSymbol, R_NilValue))));
			nProtected ++;
		}

		// ***************************************************************
		// for-loop calling

		bool ifend = false;
		int ans_index = 0;
		do {
			switch (VarIdx)
			{
				case 2:
					INTEGER(R_Index)[0] = ans_index + 1;
					break;
				case 3:
					INTEGER(R_Index)[0] = NodeList.begin()->_Index + 1;
					break;
			}
			if (NodeList.size() <= 1)
			{
				// ToDo: optimize this
				SEXP tmp = NodeList[0].NeedRData(nProtected);
				if (tmp != R_call_param)
				{
					R_call_param = tmp;
					if (VarIdx > 1)
					{
						PROTECT(R_fcall = LCONS(FUN, LCONS(R_Index,
							LCONS(R_call_param, LCONS(R_DotsSymbol, R_NilValue)))));
					} else {
						PROTECT(R_fcall = LCONS(FUN,
							LCONS(R_call_param, LCONS(R_DotsSymbol, R_NilValue))));
					}
					nProtected ++;
				}
				NodeList[0].ReadData(R_call_param);
			} else {
				int idx = 0;
				for (it=NodeList.begin(); it != NodeList.end(); it ++)
				{
					SEXP tmp = it->NeedRData(nProtected);
					it->ReadData(tmp);
					SET_ELEMENT(R_call_param, idx, tmp);
					idx ++;
				}
			}

			// call R function
			SEXP val = eval(R_fcall, rho);
			switch (DatType)
			{
			case 0:    // integer
				val = AS_INTEGER(val);
				INTEGER(rv_ans)[ans_index] = (LENGTH(val) > 0) ?
					INTEGER(val)[0] : NA_INTEGER;
				break;
			case 1:    // double
				val = AS_NUMERIC(val);
				REAL(rv_ans)[ans_index] = (LENGTH(val) > 0) ?
					REAL(val)[0] : R_NaN;
				break;
			case 2:    // character
				val = AS_CHARACTER(val);
				SET_STRING_ELT(rv_ans, ans_index,
					(LENGTH(val) > 0) ? STRING_ELT(val, 0) : NA_STRING);
				break;
			case 3:    // others
				if (NAMED(val) > 0)
				{
					// the object is bound to other symbol(s), need a copy
					val = duplicate(val);
				}
				SET_ELEMENT(rv_ans, ans_index, val);
				break;
			}
			ans_index ++;

			// check the end
			for (it=NodeList.begin(); it != NodeList.end(); it ++)
			{
				if (!it->NextCell())
					{ ifend = true; break; }
			}

		} while (!ifend);

		// finally
		UNPROTECT(nProtected);

	CORECATCH_CALL

	// output
	return(rv_ans);
}



// ###########################################################
// Apply functions via a sliding window over variants
// ###########################################################

/// Apply functions via a sliding window over variants
DLLEXPORT SEXP seq_SlidingWindow(SEXP gdsfile, SEXP var_name,
	SEXP win_size, SEXP shift_size, SEXP FUN, SEXP as_is, SEXP var_index,
	SEXP rho)
{
	SEXP rv_ans = R_NilValue;
	CORETRY_CALL

		// the selection
		TInitObject::TSelection &Sel = Init.Selection(gdsfile);
		// the GDS root node
		PdGDSObj Root = GDS_OBJECT(getListElement(gdsfile, "root"));

		// init selection
		if (Sel.Sample.empty())
		{
			PdSequenceX N = CHECK(gds_NodePath(Root, "sample.id"));
			int Cnt = gds_SeqGetCount(N);
			if (Cnt < 0) throw ErrSeqArray("Invalid dimension of 'sample.id'.");
			Sel.Sample.resize(Cnt, TRUE);
		}
		if (Sel.Variant.empty())
		{
			PdSequenceX N = CHECK(gds_NodePath(Root, "variant.id"));
			int Cnt = gds_SeqGetCount(N);
			if (Cnt < 0) throw ErrSeqArray("Invalid dimension of 'variant.id'.");
			Sel.Variant.resize(Cnt, TRUE);
		}

		// the number of calling PROTECT
		int nProtected = 0;
		// the number of selected variants
		int nVariant = 0;
		for (vector<CBOOL>::iterator it = Sel.Variant.begin();
			it != Sel.Variant.end(); it ++)
		{
			if (*it) nVariant ++;
		}
		if (nVariant <= 0)
			throw ErrSeqArray("There is no selected variant.");

		// sliding window size
		int wsize = INTEGER(win_size)[0];
		if ((wsize > nVariant) || (wsize <= 0))
			throw ErrSeqArray("`win.size' is out of range (1..%d).", nVariant);

		// shift
		int shift = INTEGER(shift_size)[0];
		if (shift <= 0)
			throw ErrSeqArray("`shift' should be greater than 0.");

		// ***************************************************************
		// initialize the GDS Node list

		vector<TVariable_ApplyByVariant> NodeList(Rf_length(var_name));
		vector<TVariable_ApplyByVariant>::iterator it;

		// for - loop
		for (int i=0; i < Rf_length(var_name); i++)
		{
			// the path of GDS variable
			string s = CHAR(STRING_ELT(var_name, i));
			TVariable_ApplyByVariant::TType VarType;

			if ( s=="variant.id" || s=="position" || s=="chromosome" ||
				s=="allele" || s=="annotation/id" || s=="annotation/qual" ||
				s=="annotation/filter" )
			{
				// ***********************************************************
				// variant.id, position, chromosome, allele, annotation/id
				// annotation/qual, annotation/filter
				VarType = TVariable_ApplyByVariant::ctBasic;
			} else if (s == "genotype")
			{
				VarType = TVariable_ApplyByVariant::ctGenotype;
				s.append("/data");
			} else if (s == "phase")
			{
				// *******************************************************
				// phase/
				VarType = TVariable_ApplyByVariant::ctPhase;
				s.append("/data");
			} else if (strncmp(s.c_str(), "annotation/info/", 16) == 0)
			{
				VarType = TVariable_ApplyByVariant::ctInfo;
			} else if (strncmp(s.c_str(), "annotation/format/", 18) == 0)
			{
				VarType = TVariable_ApplyByVariant::ctFormat;
				s.append("/data");
			} else {
				throw ErrSeqArray("'%s' is not a standard variable name, and the standard format:\n"
					"\tvariant.id, position, chromosome, allele, "
					"annotation/id, annotation/qual, annotation/filter\n"
					"\tannotation/info/VARIABLE_NAME', annotation/format/VARIABLE_NAME",
					s.c_str());
			}

			NodeList[i].InitObject(VarType, s.c_str(), Root, Sel.Variant.size(),
				&Sel.Variant[0], Sel.Sample.size(), &Sel.Sample[0]);
		}

		// ***********************************************************
		// as.is
		//     0: integer, 1: double, 2: character, 3: list, other: NULL
		int DatType;
		const char *as = CHAR(STRING_ELT(as_is, 0));
		if (strcmp(as, "integer") == 0)
			DatType = 0;
		else if (strcmp(as, "double") == 0)
			DatType = 1;
		else if (strcmp(as, "character") == 0)
			DatType = 2;
		else if (strcmp(as, "list") == 0)
			DatType = 3;
		else if (strcmp(as, "none") == 0)
			DatType = -1;
		else
			throw ErrSeqArray("'as.is' is not valid!");

		// initialize the return value
		R_xlen_t new_len = (nVariant - wsize + 1);
		new_len = (new_len / shift) + ((new_len % shift) ? 1 : 0);

		switch (DatType)
		{
		case 0:
			PROTECT(rv_ans = NEW_INTEGER(new_len));
			nProtected ++;
			break;
		case 1:
			PROTECT(rv_ans = NEW_NUMERIC(new_len));
			nProtected ++;
			break;
		case 2:
			PROTECT(rv_ans = NEW_CHARACTER(new_len));
			nProtected ++;
			break;
		case 3:
			PROTECT(rv_ans = NEW_LIST(new_len));
			nProtected ++;
			break;
		default:
			rv_ans = R_NilValue;
		}

		// ***********************************************************
		// rho
		if (!isEnvironment(rho))
			throw ErrSeqArray("'rho' should be an environment");


		// ***************************************************************
		// initialize calling

		// 1 -- none, 2 -- relative, 3 -- absolute
		int VarIdx = INTEGER(var_index)[0];

		SEXP R_fcall, R_call_param, R_Index=NULL;
		PROTECT(R_call_param = NEW_LIST(wsize));
		nProtected ++;
		if (VarIdx > 1)
		{
			PROTECT(R_Index = NEW_INTEGER(1));
			nProtected ++;
			PROTECT(R_fcall = LCONS(FUN, LCONS(R_Index,
				LCONS(R_call_param, LCONS(R_DotsSymbol, R_NilValue)))));
			nProtected ++;
		} else {
			PROTECT(R_fcall = LCONS(FUN,
				LCONS(R_call_param, LCONS(R_DotsSymbol, R_NilValue))));
			nProtected ++;
		}


		// ***************************************************************
		// for-loop calling

		// initialize the sliding window
		for (int i=1; i < wsize; i++)
		{
			if (NodeList.size() > 1)
			{
				SEXP _param = NEW_LIST(NodeList.size());
				SET_ELEMENT(R_call_param, i, _param);
				SET_NAMES(_param, GET_NAMES(var_name));

				int idx = 0;
				for (it=NodeList.begin(); it != NodeList.end(); it ++)
				{
					SEXP tmp = it->NeedRData(nProtected);
					it->ReadData(tmp);
					SET_ELEMENT(_param, idx, duplicate(tmp));
					idx ++;
				}
			} else {
				SEXP tmp = NodeList[0].NeedRData(nProtected);
				NodeList[0].ReadData(tmp);
				SET_ELEMENT(R_call_param, i, duplicate(tmp));
			}

			// check the end
			for (it=NodeList.begin(); it != NodeList.end(); it ++)
			{
				if (!it->NextCell())
					throw ErrSeqArray("internal error in 'seq_SlidingWindow'");
			}
		}

		bool ifend = false;
		int ans_index, variant_index, shift_step;
		ans_index = variant_index = shift_step = 0;

		do {
			// push
			for (int i=1; i < wsize; i++)
			{
				SET_ELEMENT(R_call_param, i-1,
					VECTOR_ELT(R_call_param, i));
			}
			SET_ELEMENT(R_call_param, wsize-1, R_NilValue);

			if (NodeList.size() > 1)
			{
				SEXP _param = NEW_LIST(NodeList.size());
				SET_ELEMENT(R_call_param, wsize-1, _param);
				SET_NAMES(_param, GET_NAMES(var_name));

				int idx = 0;
				for (it=NodeList.begin(); it != NodeList.end(); it ++)
				{
					SEXP tmp = it->NeedRData(nProtected);
					it->ReadData(tmp);
					SET_ELEMENT(_param, idx, duplicate(tmp));
					idx ++;
				}
			} else {
				SEXP tmp = NodeList[0].NeedRData(nProtected);
				NodeList[0].ReadData(tmp);
				SET_ELEMENT(R_call_param, wsize-1, duplicate(tmp));
			}

			variant_index ++;


			if (shift_step <= 0)
			{
				switch (VarIdx)
				{
					case 2:
						INTEGER(R_Index)[0] = variant_index;
						break;
					case 3:
						INTEGER(R_Index)[0] = NodeList.begin()->_Index - wsize + 2;
						break;
				}

				// call R function
				SEXP val = eval(R_fcall, rho);
				switch (DatType)
				{
				case 0:    // integer
					val = AS_INTEGER(val);
					INTEGER(rv_ans)[ans_index] = (LENGTH(val) > 0) ?
						INTEGER(val)[0] : NA_INTEGER;
					break;
				case 1:    // double
					val = AS_NUMERIC(val);
					REAL(rv_ans)[ans_index] = (LENGTH(val) > 0) ?
						REAL(val)[0] : R_NaN;
					break;
				case 2:    // character
					val = AS_CHARACTER(val);
					SET_STRING_ELT(rv_ans, ans_index,
						(LENGTH(val) > 0) ? STRING_ELT(val, 0) : NA_STRING);
					break;
				case 3:    // others
					if (NAMED(val) > 0)
					{
						// the object is bound to other symbol(s), need a copy
						val = duplicate(val);
					}
					SET_ELEMENT(rv_ans, ans_index, val);
					break;
				}
				ans_index ++;
				shift_step = shift;
			}
			shift_step --;

			// check the end
			for (it=NodeList.begin(); it != NodeList.end(); it ++)
			{
				if (!it->NextCell())
					{ ifend = true; break; }
			}

		} while (!ifend);

		// finally
		UNPROTECT(nProtected);

	CORECATCH_CALL

	// output
	return(rv_ans);
}






// ###########################################################
// VCF4 parse
// ###########################################################

/// the number of alleles per site
DLLEXPORT SEXP seq_NumOfAllele(SEXP allele_node)
{
	SEXP rv_ans = R_NilValue;
	bool has_error = false;

	CORETRY

		// GDS nodes
		PdSequenceX N;
		memcpy(&N, INTEGER(allele_node), sizeof(N));

		if (gds_SeqDimCnt(N) != 1)
			throw ErrSeqArray("Invalid dimension!");
		int Count = 0;
		gds_SeqGetDim(N, &Count);

		// allocate integers
		PROTECT(rv_ans = NEW_INTEGER(Count));

		int *base = INTEGER(rv_ans);
		string s;

		for (int i=0; i < Count; i ++)
		{
			CoreArray::Int32 _st = i;
			CoreArray::Int32 _cnt = 1;
			CHECK(gds_rData(N, &_st, &_cnt, &s, svStrUTF8));

			// determine how many alleles
			int num_allele = 0;
			const char *p = s.c_str();
			while (*p != 0)
			{
				num_allele ++;
				while ((*p != 0) && (*p != ',')) p ++;
				if (*p == ',') p ++;
			}
			base[i] = num_allele;
		}

		UNPROTECT(1);

	CORECATCH(has_error = true);
	if (has_error)
		error(gds_LastError().c_str());

	// output
	return(rv_ans);
}




// ###########################################################
// analysis
// ###########################################################

DLLEXPORT SEXP seq_Merge_Pos(SEXP opfile, SEXP outgds_root)
{
/*
	// GDS nodes
	PdSequenceX Root;
	memcpy(&Root, INTEGER(outgds_root), sizeof(Root));
	PdSequenceX varPos = gds_NodePath(Root, "position");

	// for - loop
	for (int i=0; i < Rf_length(opfile); i++)
	{
		// get variable
		PdSequenceX _R_;
		memcpy(&_R_, INTEGER(VECTOR_ELT(opfile, i)), sizeof(_R_));
		PdSequenceX sPos = gds_NodePath(_R_, "position");

		// read and write
		
	}
*/
	return R_NilValue;
}



DLLEXPORT SEXP seq_missing_snp(SEXP geno)
{
	SEXP dim = getAttrib(geno, R_DimSymbol);
	int num_ploidy = INTEGER(dim)[0];
	int num_sample = INTEGER(dim)[1];
	int miss_cnt = 0;

	int *p = INTEGER(geno);
	for (int i=0; i < num_sample; i++)
	{
		int *pp = p;
		for (int j=0; j < num_ploidy; j++, pp++)
		{
			if (UInt32(*pp) > 2)
				{ miss_cnt ++; break; }
		}
		p += num_ploidy;
	}

	SEXP rv;
	PROTECT(rv = NEW_NUMERIC(1));
	REAL(rv)[0] = (double)miss_cnt / num_sample;
	UNPROTECT(1);

	return rv;
}


DLLEXPORT SEXP seq_missing_samp(SEXP geno, SEXP miss_cnt)
{
	SEXP dim = getAttrib(geno, R_DimSymbol);
	int num_ploidy = INTEGER(dim)[0];
	int num_sample = INTEGER(dim)[1];
	int *miss = INTEGER(miss_cnt);

	int *p = INTEGER(geno);
	for (int i=0; i < num_sample; i++)
	{
		int *pp = p;
		for (int j=0; j < num_ploidy; j++, pp++)
		{
			if (UInt32(*pp) > 2)
				{ miss[i] ++; break; }
		}
		p += num_ploidy;
	}

	return R_NilValue;
}


DLLEXPORT SEXP seq_allele_freq(SEXP geno)
{
	SEXP dim = getAttrib(geno, R_DimSymbol);
	int num_ploidy = INTEGER(dim)[0];
	int num_sample = INTEGER(dim)[1];
	int ref_cnt=0, valid_cnt=0;

	int *p = INTEGER(geno);
	for (int i=0; i < num_sample; i++)
	{
		int *pp = p;
		for (int j=0; j < num_ploidy; j++, pp++)
		{
			if (UInt32(*pp) == 0) ref_cnt ++;
			if (UInt32(*pp) <= 2) valid_cnt ++;
		}
		p += num_ploidy;
	}

	SEXP rv;
	PROTECT(rv = NEW_NUMERIC(1));
	REAL(rv)[0] = (double)ref_cnt / valid_cnt;
	UNPROTECT(1);

	return rv;
}





// ###########################################################
// Convert to VCF4: GDS -> VCF4
// ###########################################################

// double quote the text if it is needed
static string QuoteText(const char *p)
{
	string rv;

	rv.clear();
	bool flag = false;
	for (; *p != 0; p ++)
	{
		switch (*p)
		{
			case ',': case ';':
				flag = true; rv.push_back(*p); break;
			case '\"':
				flag = true; rv.append("\\\""); break;
			case '\'':
				flag = true; rv.append("\\\'"); break;
			case ' ':
				flag = true; rv.push_back(' '); break;
			default:
				rv.push_back(*p);
		}
	}
	if (flag) // add double quote
	{
		rv.insert(0, "\"");
		rv.push_back('\"');
	}

	return rv;
}

/// double quote text if needed
DLLEXPORT SEXP seq_Quote(SEXP text, SEXP dQuote)
{
	SEXP NewText, ans;
	PROTECT(NewText = AS_CHARACTER(text));
	PROTECT(ans = NEW_CHARACTER(Rf_length(NewText)));

	for (int i=0; i < Rf_length(NewText); i++)
	{
		string tmp = QuoteText(CHAR(STRING_ELT(NewText, i)));
		if (LOGICAL(dQuote)[0] == TRUE)
		{
			if ((tmp[0] != '\"') || (tmp[tmp.size()-1] != '\"'))
			{
				tmp.insert(0, "\"");
				tmp.push_back('\"');
			}
		}
		SET_STRING_ELT(ans, i, mkChar(tmp.c_str()));
	}

	UNPROTECT(2);
	return ans;
}


// convert to a string
static const string TO_TEXT(SEXP X, int Start=0, int MaxCnt=-1,
	bool VarLength=false, bool NoBlank=true, int Step=1)
{
	char buffer[64];
	string ans;

	if (MaxCnt < 0)
		MaxCnt = (Rf_length(X) - Start) / Step;

	if (IS_INTEGER(X) || IS_LOGICAL(X))
	{
		int *Base = (IS_INTEGER(X) ? INTEGER(X) : LOGICAL(X)) + Start;
		if (VarLength || !NoBlank)
		{
			for (; MaxCnt > 0; MaxCnt --)
				if (Base[(MaxCnt-1)*Step] != NA_INTEGER) break;
		}
		for (int i=0; i < MaxCnt; i++, Base += Step)
		{
			if (i > 0) ans.push_back(',');
			if (*Base != NA_INTEGER)
			{
				snprintf(buffer, sizeof(buffer), "%d", *Base);
				ans.append(buffer);
			} else
				ans.push_back('.');
		}
	} else if (IS_NUMERIC(X))
	{
		double *Base = REAL(X) + Start;
		if (VarLength || !NoBlank)
		{
			for (; MaxCnt > 0; MaxCnt --)
				if (R_finite(Base[(MaxCnt-1)*Step])) break;
		}
		for (int i=0; i < MaxCnt; i++, Base += Step)
		{
			if (i > 0) ans.push_back(',');
			if (R_finite(*Base))
			{
				snprintf(buffer, sizeof(buffer), "%0.6g", *Base);
				ans.append(buffer);
			} else
				ans.push_back('.');
		}
	} else if (IS_CHARACTER(X) || Rf_isFactor(X))
	{
		if (Rf_isFactor(X))
			X = Rf_asCharacterFactor(X);
		if (VarLength || !NoBlank)
		{
			for (; MaxCnt > 0; MaxCnt --)
			{
				SEXP s = STRING_ELT(X, Start + (MaxCnt-1)*Step);
				if ((s != NA_STRING) && (CHAR(s)[0] != 0)) break;
			}
		}
		for (int i=0; i < MaxCnt; i ++, Start += Step)
		{
			if (i > 0) ans.push_back(',');
			if (STRING_ELT(X, Start) != NA_STRING)
				ans.append(QuoteText(CHAR(STRING_ELT(X, Start))));
			else
				ans.push_back('.');
		}
	}

	if (NoBlank)
	{
		if (ans.empty()) ans = ".";
	}

	return ans;
}


/// used in seq_OutVCF4
static vector<int> _VCF4_INFO_Number;    //< 
static vector<int> _VCF4_FORMAT_Number;  //< 

/// convert to VCF4
DLLEXPORT SEXP seq_InitOutVCF4(SEXP Info, SEXP Format)
{
	int *pInfo = INTEGER(Info);
	_VCF4_INFO_Number.assign(pInfo, pInfo + Rf_length(Info));
	int *pFmt = INTEGER(Format);
	_VCF4_FORMAT_Number.assign(pFmt, pFmt + Rf_length(Format));

	return R_NilValue;
}

/// convert to VCF4
DLLEXPORT SEXP seq_OutVCF4(SEXP X)
{
	const char *p, *s;
	string txt, tmp;
	int n;

	// variable list
	SEXP VarNames = getAttrib(X, R_NamesSymbol);

	// *************************************************************************
	// the first seven columns: chr, pos, id, allele (REF/ALT), qual, filter

	// CHROM
	txt.append(TO_TEXT(VECTOR_ELT(X, 0)));
	txt.push_back('\t');
	// POS
	txt.append(TO_TEXT(VECTOR_ELT(X, 1)));
	txt.push_back('\t');
	// ID
	txt.append(TO_TEXT(VECTOR_ELT(X, 2)));
	txt.push_back('\t');

	// allele -- REF/ALT
	s = p = CHAR(STRING_ELT(AS_CHARACTER(VECTOR_ELT(X, 3)), 0));
	n = 0;
	while ((*p != 0) && (*p != ','))
	{
		n ++; p ++;
	}

	// REF
	if (n > 0) txt.append(s, n); else txt.push_back('.');
	txt.push_back('\t');
	// ALT
	if (*p != 0)
	{
		p ++; txt.append((*p) ? p : ".");
	} else {
		txt.push_back('.');
	}
	txt.push_back('\t');

	// QUAL
	txt.append(TO_TEXT(VECTOR_ELT(X, 4)));
	txt.push_back('\t');
	// FILTER
	txt.append(TO_TEXT(VECTOR_ELT(X, 5)));
	txt.push_back('\t');


	// *************************************************************************
	// INFO

	bool NeedSeparator = false;
	n = 0;
	for (int i=0; i < (int)_VCF4_INFO_Number.size(); i++)
	{
		// name
		const char *nm = CHAR(STRING_ELT(VarNames, i + 8)) + 5;
		// SEXP
		SEXP D = VECTOR_ELT(X, i + 8);

		if (IS_LOGICAL(D))  // FLAG type
		{
			if (LOGICAL(D)[0] == TRUE)
			{
				if (NeedSeparator) txt.push_back(';');
				NeedSeparator = true;
				txt.append(nm);
				n ++;
			}
		} else {
			int L = _VCF4_INFO_Number[i];
			tmp = TO_TEXT(D, 0, (L < 0) ? -1 : L, (L < 0), false);
			if (!tmp.empty())
			{
				if (NeedSeparator) txt.push_back(';');
				NeedSeparator = true;
				txt.append(nm);
				txt.push_back('=');
				txt.append(tmp);
				n ++;
			}
		}
	}
	if (n <= 0) txt.push_back('.');	
	txt.push_back('\t');


	// *************************************************************************
	// FORMAT

	vector< pair<SEXP, int> > fmt_list;
	txt.append("GT");
	for (int i=0; i < (int)_VCF4_FORMAT_Number.size(); i ++)
	{
		const char *nm = CHAR(STRING_ELT(VarNames, i + 8 + _VCF4_INFO_Number.size()));
		SEXP D = VECTOR_ELT(X, i + 8 + _VCF4_INFO_Number.size());
		if (!isNull(D))
		{
			txt.push_back(':');
			txt.append(nm + 4);
			fmt_list.push_back(pair<SEXP, int>(D, i));
		}
	}
	txt.push_back('\t');


	// *************************************************************************
	// Genotypic data

	// genotype
	SEXP geno = VECTOR_ELT(X, 6);
	SEXP geno_dim = GET_DIM(geno);
	if (Rf_length(geno_dim) != 2)
		error("Invalid dimension of genotypes.");

	const int NumAllele = INTEGER(geno_dim)[0];
	const int NumSample = INTEGER(geno_dim)[1];

	// phase information
	SEXP phase = VECTOR_ELT(X, 7);
	if (Rf_length(phase) != (NumAllele-1)*NumSample)
		error("Invalid dimension of phasing information.");

	int *pSamp = INTEGER(geno);
	int *pAllele = INTEGER(phase);

	// for-loop of samples
	for (int i=0; i < NumSample; i ++)
	{
		// genotypes
		for (int j=0; j < NumAllele; j++, pSamp++)
		{
			if (j > 0)
			{
				txt.push_back(*pAllele ? '|' : '/');
				pAllele ++;
			}
			if (*pSamp != NA_INTEGER)
			{
				char buf[32];
				snprintf(buf, sizeof(buf), "%d", *pSamp);
				txt.append(buf);
			} else
				txt.push_back('.');
		}

		// annotation
		vector< pair<SEXP, int> >::iterator it;
		for (it=fmt_list.begin(); it != fmt_list.end(); it ++)
		{
			txt.push_back(':');
			int nTotal = Rf_length(it->first);
			int nColumn = nTotal / NumSample;
			if ((nTotal % NumSample) != 0)
				error("Internal Error: invalid dimension.");

			int L = _VCF4_FORMAT_Number[it->second];
			tmp = (L < 0) ? TO_TEXT(it->first, i, nColumn, true, true, NumSample) :
				TO_TEXT(it->first, i, L, false, true, NumSample);
			txt.append(tmp);
		}

		// add '\t'
		if (i < (NumSample-1)) txt.push_back('\t');
	}


	// append '\n'
	txt.push_back('\n');

	// return
	SEXP ans;
	PROTECT(ans = NEW_CHARACTER(1));
	SET_STRING_ELT(ans, 0, mkChar(txt.c_str()));
	UNPROTECT(1);
	return ans;
}






// ###########################################################
// Convert from VCF4: VCF4 -> GDS
// ###########################################################

/// trim blank characters
static void _Trim_(string &val)
{
	const char *st = val.c_str();
	const char *p = st;
	while ((*p == ' ') || (*p == '\t')) p ++;
	if (p != st) val.erase(0, p - st);

	int L = val.size();
	p = val.c_str() + L - 1;
	while (L > 0)
	{
		if ((*p != ' ') && (*p != '\t')) break;
		L --; p --;
	}
	val.resize(L);
}

/// get an integer from a string
static CoreArray::Int32 getInt32(const string &txt, bool RaiseError)
{
	const char *p = SKIP(txt.c_str());
	char *endptr = (char*)p;
	long int val = strtol(p, &endptr, 10);

	if (endptr == p)
	{
		if ((*p != '.') && RaiseError)
		{
			throw ErrSeqArray("Invalid integer conversion \"%s\".",
				SHORT_TEXT(p).c_str());
		}
		val = NA_INTEGER;
	} else {
		if ((val < INT_MIN) || (val > INT_MAX))
		{
			val = NA_INTEGER;
			if (RaiseError)
			{
				throw ErrSeqArray("Invalid integer conversion \"%s\".",
					SHORT_TEXT(p).c_str());
			}
		}
		p = SKIP(endptr);
		if (*p != 0)
		{
			val = NA_INTEGER;
			if (RaiseError)
			{
				throw ErrSeqArray("Invalid integer conversion \"%s\".",
					SHORT_TEXT(p).c_str());
			}
		}
	}
	return val;
}

/// get multiple integers from a string
static void getInt32Array(const string &txt, vector<CoreArray::Int32> &I32,
	bool RaiseError)
{
	const char *p = SKIP(txt.c_str());
	I32.clear();
	while (*p)
	{
		char *endptr = (char*)p;
		long int val = strtol(p, &endptr, 10);
		
		if (endptr == p)
		{
			if ((*p != '.') && RaiseError)
			{
				throw ErrSeqArray("Invalid integer conversion \"%s\".",
					SHORT_TEXT(p).c_str());
			}
			val = NA_INTEGER;
		} else {
			if ((val < INT_MIN) || (val > INT_MAX))
			{
				val = NA_INTEGER;
				if (RaiseError)
				{
					throw ErrSeqArray("Invalid integer conversion \"%s\".",
						SHORT_TEXT(p).c_str());
				}
			}
			p = endptr;
		}

		I32.push_back(val);
		while ((*p != 0) && (*p != ',')) p ++;
		if (*p == ',') p ++;
	}
}


/// get a float number from a string
static float getFloat(string &txt, bool RaiseError)
{
	const char *p = SKIP(txt.c_str());
	char *endptr = (char*)p;
	float val = strtof(p, &endptr);
	if (endptr == p)
	{
		if ((*p != '.') && RaiseError)
		{
			throw ErrSeqArray("Invalid float conversion \"%s\".",
				SHORT_TEXT(p).c_str());
		}
		val = R_NaN;
	} else {
		p = SKIP(endptr);
		if (*p != 0)
		{
			val = R_NaN;
			if (RaiseError)
			{
				throw ErrSeqArray("Invalid float conversion \"%s\".",
					SHORT_TEXT(p).c_str());
			}
		}
	}
	return val;
}

/// get an integer from  a string
static void getFloatArray(const string &txt, vector<float> &F32,
	bool RaiseError)
{
	const char *p = SKIP(txt.c_str());
	F32.clear();
	while (*p)
	{
		char *endptr = (char*)p;
		float val = strtof(p, &endptr);
		if (endptr == p)
		{
			if ((*p != '.') && RaiseError)
			{
				throw ErrSeqArray("Invalid float conversion \"%s\".",
					SHORT_TEXT(p).c_str());
			}
			val = R_NaN;
		} else
			p = endptr;

		F32.push_back(val);
		while ((*p != 0) && (*p != ',')) p ++;
		if (*p == ',') p ++;
	}
}

/// get an integer from  a string
static void getStringArray(string &txt, vector<string> &UTF8s)
{
	string val;
	const char *p = txt.c_str();
	while ((*p == ' ') || (*p == '\t')) p ++;

	UTF8s.clear();
	while (*p)
	{
		val.clear();
		while ((*p != 0) && (*p != ','))
			{ val.push_back(*p); p ++; }
		_Trim_(val);
		UTF8s.push_back(val);
		if (*p == ',') p ++;
	}
}

/// get name and value
static const char *_GetNameValue(const char *p, string &name, string &val)
{
	// name = val
	name.clear();
	while ((*p != 0) && (*p != ';') && (*p != '='))
	{
		name.push_back(*p);
		p ++;
	}

	val.clear();
	if (*p == '=') p ++;
	while ((*p != 0) && (*p != ';'))
	{
		val.push_back(*p);
		p ++;
	}

	if (*p == ';') p ++;
	return p;
}


/// VCF4 --> GDS
DLLEXPORT SEXP seq_Parse_VCF4(SEXP vcf_fn, SEXP header, SEXP gds_root,
	SEXP param, SEXP ReadLineFun, SEXP ReadLine_Param, SEXP rho)
{
	SEXP rv_ans = R_NilValue;
	bool has_error = false;
	const char *fn = CHAR(STRING_ELT(vcf_fn, 0));

	// define a variable for reading lines
	CReadLine RL;

	// cell buffer
	string cell;
	cell.reserve(4096);


	CORETRY

		// the number of calling PROTECT
		int nProtected = 0;


		// *********************************************************
		// initialize variables		

		// the total number of samples
		int nTotalSamp = INTEGER(getListElement(param, "sample.num"))[0];
		// the variable name for genotypic data
		string geno_id = CHAR(STRING_ELT(getListElement(param, "genotype.var.name"), 0));
		// NumericRaiseError
		bool RaiseError = (LOGICAL(getListElement(param, "raise.error"))[0] == TRUE);
		// verbose
		// bool Verbose = (LOGICAL(getListElement(param, "verbose"))[0] == TRUE);

		// the number of ploidy
		int num_ploidy = INTEGER(getListElement(header, "num.ploidy"))[0];
		if (num_ploidy <= 0)
			throw ErrSeqArray("Invalid header$num.ploidy: %d.", num_ploidy);


		// GDS nodes
		PdSequenceX Root;
		memcpy(&Root, INTEGER(gds_root), sizeof(Root));

		PdSequenceX varIdx = gds_NodePath(Root, "variant.id");
		PdSequenceX varChr = gds_NodePath(Root, "chromosome");
		PdSequenceX varPos = gds_NodePath(Root, "position");
		PdSequenceX varRSID = gds_NodePath(Root, "annotation/id");
		PdSequenceX varAllele = gds_NodePath(Root, "allele");

		PdSequenceX varQual = gds_NodePath(Root, "annotation/qual");
		PdSequenceX varFilter = gds_NodePath(Root, "annotation/filter");

		PdSequenceX varGeno = gds_NodePath(Root, "genotype/data");
		PdSequenceX varGenoLen = gds_NodePath(Root, "genotype/@data");
		PdSequenceX varGenoExtraIdx = gds_NodePath(Root, "genotype/extra.index");
		PdSequenceX varGenoExtra = gds_NodePath(Root, "genotype/extra");

		PdSequenceX varPhase = gds_NodePath(Root, "phase/data");
		PdSequenceX varPhaseExtraIdx = gds_NodePath(Root, "phase/extra.index");
		PdSequenceX varPhaseExtra = gds_NodePath(Root, "phase/extra");


		// getListElement: info
		vector<TVCF_Field_Info> info_list;
		{
			SEXP info = getListElement(header, "info");
			SEXP info_ID = getListElement(info, "ID");
			SEXP info_inttype = getListElement(info, "int_type");
			SEXP info_intnum = getListElement(info, "int_num");
			SEXP info_flag = getListElement(info, "import.flag");
			TVCF_Field_Info val;

			for (int i=0; i < Rf_length(info_ID); i++)
			{
				val.name = CHAR(STRING_ELT(info_ID, i));
				val.type = INTEGER(info_inttype)[i];
				val.import_flag = (LOGICAL(info_flag)[i] == TRUE);
				val.number = INTEGER(info_intnum)[i];
				val.data_obj = gds_NodePath(Root,
					(string("annotation/info/") + val.name).c_str());
				val.len_obj = gds_NodePath(Root,
					(string("annotation/info/@") + val.name).c_str());
				info_list.push_back(val);
			}
		}

		// getListElement: format
		vector<TVCF_Field_Format> format_list;
		{
			SEXP fmt = getListElement(header, "format");
			SEXP fmt_ID = getListElement(fmt, "ID");
			SEXP fmt_inttype = getListElement(fmt, "int_type");
			SEXP fmt_intnum = getListElement(fmt, "int_num");
			SEXP fmt_flag = getListElement(fmt, "import.flag");
			TVCF_Field_Format val;

			for (int i=0; i < Rf_length(fmt_ID); i++)
			{
				val.name = CHAR(STRING_ELT(fmt_ID, i));
				val.type = INTEGER(fmt_inttype)[i];
				val.import_flag = (LOGICAL(fmt_flag)[i] == TRUE);
				val.number = INTEGER(fmt_intnum)[i];
				val.data_obj = gds_NodePath(Root,
					(string("annotation/format/") + val.name + "/data").c_str());
				val.len_obj = gds_NodePath(Root,
					(string("annotation/format/") + val.name + "/@data").c_str());
				format_list.push_back(val);

				switch (val.type)
				{
					case FIELD_TYPE_INT:
						format_list.back().I32ss.resize(nTotalSamp);
						break;
					case FIELD_TYPE_FLOAT:
						format_list.back().F32ss.resize(nTotalSamp);
						break;
					case FIELD_TYPE_STRING:
						format_list.back().UTF8ss.resize(nTotalSamp);
						break;
					default:
						throw ErrSeqArray("Invalid FORMAT Type.");
				}
			}
		}

		// filter string list
		vector<string> filter_list;

		// variant id (integer)
		CoreArray::Int32 variant_index = 1;

		// the string buffer
		string name, value;
		name.reserve(256);
		value.reserve(4096);

		// the numeric buffer
		CoreArray::Int32 I32;
		vector<CoreArray::Int32> I32s;
		I32s.reserve(nTotalSamp);

		CoreArray::Float32 F32;
		vector<CoreArray::Float32> F32s;
		F32s.reserve(nTotalSamp);

		// the string buffer
		vector<string> StrList;
		StrList.reserve(nTotalSamp);

		// genotypes
		vector< vector<Int16> > Geno;
		Geno.resize(nTotalSamp);
		vector<Int8> I8s;
		I8s.reserve(nTotalSamp * num_ploidy);

		const char *pCh;
		vector<TVCF_Field_Info>::iterator pInfo;
		vector< TVCF_Field_Format* >::iterator pFormat;
		vector< TVCF_Field_Format* > fmt_ptr;
		fmt_ptr.reserve(format_list.size());

		// the number of alleles in total at a specified site
		int num_allele;


		// *********************************************************
		// initialize calling
		SEXP R_Read_Call;
		PROTECT(R_Read_Call = LCONS(ReadLineFun, LCONS(ReadLine_Param, R_NilValue)));
		nProtected ++;
		RL.Init(R_Read_Call, rho);


		// *********************************************************
		// skip the header

		while (!RL.IfEnd())
		{
			const char *p = RL.ReadLine();
			if (strncmp(p, "#CHROM", 6) == 0)
				break;
		}


		// *********************************************************
		// parse the context

		while (!RL.IfEnd())
		{
			// *****************************************************
			// scan line by line

			// #####################################################
			// variant id
			CHECK(gds_AppendData(varIdx, 1, &variant_index, svInt32));
			variant_index ++;


			// #####################################################
			// column 1: CHROM
			RL.GetCell(cell, false);
			CHECK(gds_AppendString(varChr, cell.c_str()));


			// #####################################################
			// column 2: POS
			RL.GetCell(cell, false);
			I32 = getInt32(cell, RaiseError);
			CHECK(gds_AppendData(varPos, 1, &I32, svInt32));


			// #####################################################
			// column 3: ID
			RL.GetCell(cell, false);
			if (cell == ".") cell.clear();
			CHECK(gds_AppendString(varRSID, cell.c_str()));


			// #####################################################
			// column 4 & 5: REF + ALT 
			RL.GetCell(value, false);
			RL.GetCell(cell, false);
			if (!cell.empty() && (cell != "."))
			{
				value.push_back(',');
				value.append(cell);
			}
			CHECK(gds_AppendString(varAllele, value.c_str()));
			// determine how many alleles
			num_allele = 0;
			pCh = value.c_str();
			while (*pCh != 0)
			{
				num_allele ++;
				while ((*pCh != 0) && (*pCh != ',')) pCh ++;
				if (*pCh == ',') pCh ++;
			}


			// #####################################################
			// column 6: QUAL
			RL.GetCell(cell, false);
			F32 = getFloat(cell, RaiseError);
			CHECK(gds_AppendData(varQual, 1, &F32, svFloat32));


			// #####################################################
			// column 7: FILTER
			RL.GetCell(cell, false);
			if (!cell.empty() && (cell != "."))
			{
				vector<string>::iterator p =
					find(filter_list.begin(), filter_list.end(), cell);
				if (p == filter_list.end())
				{
					filter_list.push_back(cell);
					I32 = filter_list.size();
				} else
					I32 = p - filter_list.begin() + 1;
			} else
				I32 = NA_INTEGER;
			CHECK(gds_AppendData(varFilter, 1, &I32, svInt32));


			// #####################################################
			// column 8: INFO

			// initialize
			for (pInfo = info_list.begin(); pInfo != info_list.end(); pInfo++)
				pInfo->used = false;
			RL.GetCell(cell, false);

			// parse
			pCh = cell.c_str();
			while (*pCh)
			{
				pCh = _GetNameValue(pCh, name, value);
				for (pInfo=info_list.begin(); pInfo != info_list.end(); pInfo++)
				{
					if (pInfo->name == name)
						break;
				}
				if (pInfo != info_list.end())
				{
					if (pInfo->used)
						throw ErrSeqArray("Duplicate INFO ID: %s.", name.c_str());

					if (pInfo->import_flag)
					{
						switch (pInfo->type)
						{
						case FIELD_TYPE_INT:
							getInt32Array(value, I32s, RaiseError);
							pInfo->Check(I32s, name, num_allele);
							CHECK(gds_AppendData(pInfo->data_obj, I32s.size(),
								&(I32s[0]), svInt32));
							break;

						case FIELD_TYPE_FLOAT:
							getFloatArray(value, F32s, RaiseError);
							pInfo->Check(F32s, name, num_allele);
							CHECK(gds_AppendData(pInfo->data_obj, F32s.size(),
								&(F32s[0]), svFloat32));
							break;

						case FIELD_TYPE_FLAG:
							if (!value.empty())
							{
								throw ErrSeqArray("INFO ID '%s' should be a flag without values.",
									name.c_str());
							}
							I32 = 1;
							CHECK(gds_AppendData(pInfo->data_obj, 1, &I32, svInt32));
							break;

						case FIELD_TYPE_STRING:
							getStringArray(value, StrList);
							pInfo->Check(StrList, name, num_allele);
							for (int k=0; k < (int)StrList.size(); k++)
							{
								CHECK(gds_AppendString(pInfo->data_obj,
									StrList[k].c_str()));
							}
							break;

						default:
							throw ErrSeqArray("Invalid INFO Type.");
						}
					}

					pInfo->used = true;
				} else {
					if (RaiseError)
					{
						throw ErrSeqArray(
							"Unknown INFO ID: %s, should be defined ahead.",
							name.c_str());
					} else
						warning("Unknown INFO ID '%s' is ignored.", name.c_str());
				}
			}

			// for which does not exist
			for (pInfo = info_list.begin(); pInfo != info_list.end(); pInfo++)
			{
				if (!pInfo->used && pInfo->import_flag)
				{
					switch (pInfo->type)
					{
					case FIELD_TYPE_INT:
						pInfo->Fill(I32s, NA_INTEGER);
						break;
					case FIELD_TYPE_FLOAT:
						pInfo->Fill(F32s, (float)R_NaN);
						break;
					case FIELD_TYPE_FLAG:
						I32 = 0;
						CHECK(gds_AppendData(pInfo->data_obj, 1, &I32, svInt32));
						break;
					case FIELD_TYPE_STRING:
						pInfo->Fill(StrList, string());
						break;
					default:
						throw ErrSeqArray("Invalid INFO Type.");
					}
					pInfo->used = true;
				}
			}


			// #####################################################
			// column 9: FORMAT

			// initialize
			for (pFormat = fmt_ptr.begin(); pFormat != fmt_ptr.end(); pFormat++)
				(*pFormat)->used = false;
			RL.GetCell(cell, false);

			// parse
			bool first_id_flag = true;
			fmt_ptr.clear();
			pCh = cell.c_str();
			while (*pCh)
			{
				name.clear();
				while ((*pCh != 0) && (*pCh != ':'))
					{ name.push_back(*pCh); pCh ++; }
				if (*pCh == ':') pCh ++;
				_Trim_(name);

				if (first_id_flag)
				{
					// genotype ID
					if (name != geno_id)
					{
						throw ErrSeqArray("The first FORMAT ID should be '%s'.",
							geno_id.c_str());
					}
					first_id_flag = false;

				} else {
					// find ID
					vector<TVCF_Field_Format>::iterator it;
					for (it = format_list.begin(); it != format_list.end(); it++)
					{
						if (it->name == name)
							{ it->used = true; break; }
					}
					if (it == format_list.end())
					{
						if (RaiseError)
						{
							throw ErrSeqArray(
								"Unknown FORMAT ID: %s, it should be defined ahead.",
								name.c_str());
						} else {
							warning("Unknown FORMAT ID '%s' is ignored.",
								name.c_str());
							fmt_ptr.push_back(NULL);
						}
					} else {
						// push
						fmt_ptr.push_back(&(*it));
					}
				}
			}


			// #####################################################
			// columns for samples

			// for-loop
			for (int samp_idx=0; samp_idx < nTotalSamp; samp_idx ++)
			{
				const char *p;

				// read
				RL.GetCell(cell, samp_idx >= (nTotalSamp-1));

				// #################################################
				// the first field -- genotypes
				pCh = p = cell.c_str();
				while ((*p != 0) && (*p != ':')) p ++;
				value.assign(pCh, p);
				pCh = (*p == ':') ? (p + 1) : p;

				I32s.clear();
				vector<Int16> &pAllele = Geno[samp_idx];
				pAllele.clear();

				p = SKIP(value.c_str());
				while (*p)
				{
					char *endptr = (char*)p;
					CoreArray::Int32 val = strtol(p, &endptr, 10);

					if (endptr == p)
					{
						if ((*p != '.') && RaiseError)
							throw ErrSeqArray("Invalid integer conversion \"%s\".", SHORT_TEXT(p).c_str());
						val = -1;
					} else {
						if (val < 0)
						{
							val = -1;
							if (RaiseError)
								throw ErrSeqArray("Genotype code should be non-negative \"%s\".", SHORT_TEXT(p).c_str());
						} else if (val >= num_allele)
						{
							val = -1;
							if (RaiseError)
								throw ErrSeqArray("Genotype code is out of range \"%s\".", SHORT_TEXT(p).c_str());
						}
						p = endptr;
					}

					pAllele.push_back(val);
					while ((*p != 0) && (*p != '|') && (*p != '/'))
						p ++;
					if (*p == '|')
					{
						I32s.push_back(1); p ++;
					} else if (*p == '/')
					{
						I32s.push_back(0); p ++;
					}
				}

				// check pAllele
				if ((int)pAllele.size() < num_ploidy)
					pAllele.resize(num_ploidy, -1);

				// write phasing information
				if ((int)I32s.size() < (num_ploidy-1))
					I32s.resize(num_ploidy-1, 0);
				CHECK(gds_AppendData(varPhase, num_ploidy-1, &(I32s[0]), svInt32));
				if ((int)I32s.size() > (num_ploidy-1))
				{
					// E.g., triploid call: 0/0/1
					int Len = num_ploidy - (int)I32s.size() - 1;
					CHECK(gds_AppendData(varPhaseExtra, Len, &(I32s[num_ploidy-1]), svInt32));
					I32 = samp_idx + 1;
					CHECK(gds_AppendData(varPhaseExtraIdx, 1, &I32, svInt32));
					I32 = variant_index;
					CHECK(gds_AppendData(varPhaseExtraIdx, 1, &I32, svInt32));
					I32 = Len;
					CHECK(gds_AppendData(varPhaseExtraIdx, 1, &I32, svInt32));
				}

				// #################################################
				// the other field -- format id
				for (int i=0; i < (int)fmt_ptr.size(); i++)
				{
					TVCF_Field_Format *pFmt = fmt_ptr[i];
					p = pCh;
					while ((*p != 0) && (*p != ':')) p ++;

					if ((pFmt!=NULL) && pFmt->import_flag)
					{
						// get the field context
						value.assign(pCh, p);

						// parse the field
						switch (pFmt->type)
						{
						case FIELD_TYPE_INT:
							getInt32Array(value, pFmt->I32ss[samp_idx], RaiseError);
							pFmt->Check(pFmt->I32ss[samp_idx], name, num_allele, NA_INTEGER);
							break;

						case FIELD_TYPE_FLOAT:
							getFloatArray(value, pFmt->F32ss[samp_idx], RaiseError);
							pFmt->Check(pFmt->F32ss[samp_idx], name, num_allele, (float)R_NaN);
							break;

						case FIELD_TYPE_STRING:
							getStringArray(value, pFmt->UTF8ss[samp_idx]);
							pFmt->Check(pFmt->UTF8ss[samp_idx], name, num_allele, BlackString);
							break;

						default:
							throw ErrSeqArray("Invalid FORMAT Type.");
						}
					}

					pCh = (*p == ':') ? (p + 1) : p;
				}
			}

			// #################################################
			// write genotypes

			// determine how many bits
			int num_bits = 2;
			// plus ONE for missing value
			while ((num_allele + 1) > (1 << num_bits))
				num_bits += 2;
			I32 = num_bits / 2;
			CHECK(gds_AppendData(varGenoLen, 1, &I32, svInt32));

			// write to the variable "genotype"
			for (int bits=0; bits < num_bits; bits += 2)
			{
				I8s.clear();
				for (int i=0; i < nTotalSamp; i++)
				{
					vector<Int16> &pAllele = Geno[i];
					for (int j=0; j < num_ploidy; j++)
						I8s.push_back((pAllele[j] >> bits) & 0x03);
				}
				CHECK(gds_AppendData(varGeno, nTotalSamp*num_ploidy, &(I8s[0]), svInt8));
			}

			// write to "genotype/extra"
			for (int i=0; i < nTotalSamp; i++)
			{
				vector<Int16> &pGeno = Geno[i];
				if ((int)pGeno.size() > num_ploidy)
				{
					// E.g., triploid call: 0/0/1
					int Len = num_ploidy - (int)pGeno.size() - 1;
					CHECK(gds_AppendData(varGenoExtra, Len, &(pGeno[num_ploidy-1]), svInt32));
					I32 = i + 1;
					CHECK(gds_AppendData(varGenoExtraIdx, 1, &I32, svInt32));
					I32 = variant_index;
					CHECK(gds_AppendData(varGenoExtraIdx, 1, &I32, svInt32));
					I32 = Len;
					CHECK(gds_AppendData(varGenoExtraIdx, 1, &I32, svInt32));
				}
			}


			// #################################################
			// for-loop all format IDs: write
			for (vector<TVCF_Field_Format>::iterator it = format_list.begin();
				it != format_list.end(); it++)
			{
				if (it->import_flag)
				{
					if (it->used)
					{
						if (it->number > 0)
						{
							// fixed-length array
							it->WriteFixedLength();
							I32 = 1;
						} else if (it->number < 0)
						{
							// variable-length array
							I32 = it->WriteVariableLength(nTotalSamp, I32s, F32s);
						} else
							throw ErrSeqArray("Invalid FORMAT Number.");
						CHECK(gds_AppendData(it->len_obj, 1, &I32, svInt32));
					} else {
						I32 = 0;
						CHECK(gds_AppendData(it->len_obj, 1, &I32, svInt32));
					}
				}
			}
		}

		// set returned value: levels(filter)
		PROTECT(rv_ans = NEW_CHARACTER(filter_list.size()));
		for (int i=0; i < (int)filter_list.size(); i++)
			SET_STRING_ELT(rv_ans, i, mkChar(filter_list[i].c_str()));
		nProtected ++;


		UNPROTECT(nProtected);

	CORECATCH({
		char buf[4096];
		if (RL.ColumnNo() > 0)
		{
			snprintf(buf, sizeof(buf), "\tFILE: %s\n\tLINE: %d, COLUMN: %d, %s\n\t%s",
				fn, RL.LineNo(), RL.ColumnNo(), cell.c_str(), gds_LastError().c_str());
		} else {
			snprintf(buf, sizeof(buf), "\tFILE: %s\n\tLINE: %d\n\t%s",
				fn, RL.LineNo(), gds_LastError().c_str());
		}
		gds_LastError() = buf;
		has_error = true;
	});
	if (has_error)
		error(gds_LastError().c_str());

	// output
	return(rv_ans);
}


} // extern "C"
