// ===========================================================
//
// ReadBySample.h: Read data sample by sample
//
// Copyright (C) 2015    Xiuwen Zheng
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

#include "Common.h"


/// 
class COREARRAY_DLL_LOCAL CVarApplyBySample: public CVarApply
{
protected:
	PdSequenceX Node;       ///< the GDS variable

	C_Int32 VariantStart;   ///< start index according to the variants
	C_Int32 VariantCount;   ///< the length according to the variants
	size_t CellCount;       ///< the number of entries for the current sample
	vector<C_UInt8> VariantCellCnt;  ///< 
	map<size_t, SEXP> VarList;       ///< a list of SEXP variables

	C_SVType SVType;        ///< data type for GDS reading
	C_BOOL *SelPtr[3];      ///< pointers to selection
	C_BOOL *SampleSelect;   ///< pointer to sample selection

	vector<C_BOOL> VariantSelectBuffer;  ///< the buffer of selection

public:
	TType VarType;          ///< VCF data type
	int TotalNum_Sample;    ///< the total number of samples
	int Num_Variant;        ///< the number of selected variants

	int DimCnt;             ///< the number of dimensions
	C_Int32 DLen[4];        ///< the dimension size
	C_Int32 CurIndex;       ///< the index of variant, starting from ZERO


	CVarApplyBySample();

	void InitObject(TType Type, const char *Path, PdGDSObj Root,
		int nVariant, C_BOOL *VariantSel, int nSample, C_BOOL *SampleSel);
	void ResetObject();

	bool NextCell();

	void ReadGenoData(int *Base);

	void ReadData(SEXP Val);

	SEXP NeedRData(int &nProtected);
};



extern "C"
{
/// Apply functions over margins on a working space
COREARRAY_DLL_EXPORT SEXP sqa_Apply_Sample(SEXP gdsfile, SEXP var_name,
	SEXP FUN, SEXP as_is, SEXP var_index, SEXP rho);

} // extern "C"
