#if !defined MML_CONSOLE_PRINTER_H
#define MML_CONSOLE_PRINTER_H

#include "MMLBase.h"

namespace MML
{
	class ColDesc
	{
	public:
		std::string _name;
		int _width;
		int _precision;
		char _format;

		ColDesc(std::string name) :
			_name(name), _width(12), _precision(5), _format('D') {
		}
		ColDesc(std::string name, int width, int precision, char format) :
			_name(name), _width(width), _precision(precision), _format(format) { }
	};

	class VerticalVectorPrinter
	{
		std::vector<Vector<Real>*> listVec;
		std::vector<ColDesc> listColDesc;
		int _maxVecLen;

	public:
		VerticalVectorPrinter() : _maxVecLen(0) {}
		VerticalVectorPrinter(std::vector<Vector<Real>*> vecs) : listVec(vecs)
		{
			// calculate maximum length for these vectors
			_maxVecLen = 0;
			for (size_t i = 0; i < listVec.size(); i++)
			{
				if (listVec[i]->size() > _maxVecLen)
					_maxVecLen = listVec[i]->size();
			}
		}
		VerticalVectorPrinter(std::vector<ColDesc> vecNames, std::vector<Vector<Real>*> vecValues) : listColDesc(vecNames), listVec(vecValues)
		{
			// calculate maximum length for these vectors
			_maxVecLen = 0;
			for (size_t i = 0; i < listVec.size(); i++)
			{
				if (listVec[i]->size() > _maxVecLen)
					_maxVecLen = listVec[i]->size();
			}
		}

		void Print()
		{
			// first, save current stream state
			std::ios_base::fmtflags f(std::cout.flags());

			// change formatting
			std::cout << std::fixed;

			// printout header with vector names (if available)
			if (listColDesc.size() > 0)
			{
				for (size_t i = 0; i < listColDesc.size(); i++)
					std::cout << std::setw(listColDesc[i]._width) << listColDesc[i]._name << " ";
				std::cout << std::endl;
			}

			for (int i = 0; i < _maxVecLen; i++)
			{
				for (size_t j = 0; j < listVec.size(); j++)
				{
					if (i < listVec[j]->size())
						std::cout << std::setw(listColDesc[j]._width) << std::setprecision(listColDesc[j]._precision) << (*listVec[j])[i] << " ";
					else
						// TODO - make sure that all necessary spaces are printed
						std::cout << std::setw(listColDesc[j]._width) << " " << " ";
				}
				std::cout << std::endl;
			}

			// restore stream state
			std::cout.flags(f);
		}
	};

	// ConsolePrinterComplex
	// - Value je bazna klasa - Real, Vector, bool, string
	template<class _Tag, class _Value>
	class TablePrinter
	{
	public:
		std::string _tagName;
		std::vector<_Tag> _tags;

		int _tagWidth, _tagPrec;
		std::vector<std::tuple<int, int, char>> _listFormatSpec;

		std::vector<std::string> _valueNames;
		std::vector<std::vector<_Value>> _values;

		TablePrinter(std::string tagName, std::vector<std::string> valueNames) :
			_tagName(tagName), _valueNames(valueNames)
		{
			_tagWidth = 8;
			_tagPrec = 3;
			for (size_t i = 0; i < _valueNames.size(); i++)
			{
				_listFormatSpec.push_back(std::make_tuple(11, 5, ' '));
			}
		}

		TablePrinter(std::string tagName, int tagWidth, int tagPrec, std::vector<std::string> valueNames, std::vector<std::tuple<int, int, char>> listWidthPrec) :
			_tagName(tagName), _tagWidth(tagWidth), _tagPrec(tagPrec), _valueNames(valueNames), _listFormatSpec(listWidthPrec)
		{}

		void addRow(_Tag tag, std::vector<_Value> values)
		{
			if (values.size() != _valueNames.size())
				throw std::invalid_argument("Number of values does not match number of value names");

			_tags.push_back(tag);
			_values.push_back(values);
		}

		void Print()
		{
			std::cout << std::setw(_tagWidth) << _tagName << " ";
			for (size_t i = 0; i < _valueNames.size(); i++)
			{
				std::cout << std::setw(std::get<0>(_listFormatSpec[i])) << std::setprecision(std::get<1>(_listFormatSpec[i])) << _valueNames[i] << " ";
			}
			std::cout << std::endl;

			for (size_t i = 0; i < _tags.size(); i++)
			{

				if (std::get<2>(_listFormatSpec[0]) == 'S')
					std::cout << std::scientific << std::setw(_tagWidth) << std::setprecision(_tagPrec) << _tags[i] << " ";
				else
					std::cout << std::fixed << std::setw(_tagWidth) << std::setprecision(_tagPrec) << _tags[i] << " ";
				for (size_t j = 0; j < _values[i].size(); j++)
				{
					if (std::get<2>(_listFormatSpec[j]) == 'S')
						std::cout << std::scientific << std::setw(std::get<0>(_listFormatSpec[j])) << std::setprecision(std::get<1>(_listFormatSpec[j])) << _values[i][j] << " ";
					else
						std::cout << std::fixed << std::setw(std::get<0>(_listFormatSpec[j])) << std::setprecision(std::get<1>(_listFormatSpec[j])) << _values[i][j] << " ";
				}
				std::cout << std::endl;
			}
		}

		void Print(int tagWidth, int tagPrec, std::vector<std::pair<int, int>> listWidthPrec)
		{
			std::cout << std::fixed << std::setw(tagWidth) << _tagName << " ";
			for (size_t i = 0; i < _valueNames.size(); i++)
			{
				std::cout << std::setw(listWidthPrec[i].first) << _valueNames[i];
			}
			std::cout << std::endl;
			for (size_t i = 0; i < _tags.size(); i++)
			{
				std::cout << std::setw(tagWidth) << std::setprecision(tagPrec) << _tags[i] << " ";
				for (size_t j = 0; j < _values[i].size(); j++)
				{
					std::cout << std::setw(listWidthPrec[j].first) << std::setprecision(listWidthPrec[j].second) << _values[i][j];
				}
				std::cout << std::endl;
			}
		}
	};
}

#endif