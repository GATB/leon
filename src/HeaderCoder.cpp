/*****************************************************************************
 *   Leon: reference free compression for NGS reads
 *   A tool from the GATB (Genome Assembly Tool Box)
 *   Copyright (C) 2014  INRIA
 *   Authors: G.Benoit, G.Rizk, C.Lemaitre
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#include "HeaderCoder.hpp"

#include <bitset> //////////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! delete
/*
#define PRINT_DEBUG_ENCODER
#define PRINT_DEBUG_DECODER
*/
		

//====================================================================================
// ** AbstractHeaderCoder
//====================================================================================
AbstractHeaderCoder::AbstractHeaderCoder(Leon* leon) :
_headerSizeModel(256)
{
	cerr << "AbstractHeaderCoder::AbstractHeaderCoder - test seg 0" << endl;
	_leon = leon;
	_prevHeader = "";
	_currentHeader = "";
	/*
	for(int i=0; i<MAX_FIELD_COUNT; i++){
		_typeModel.push_back(Order0Model(HEADER_TYPE_COUNT+1));
		_fieldIndexModel.push_back(Order0Model(MAX_FIELD_COUNT));
		_fieldColumnModel.push_back(Order0Model(256));
		_misSizeModel.push_back(Order0Model(256));
		_asciiModel.push_back(Order0Model(128));
		_zeroModel.push_back(Order0Model(256));
		
		_numericSizeModel.push_back(Order0Model(10));
		_numericModels.push_back(vector<Order0Model>());
		for(int j=0; j<8; j++)
			_numericModels[i].push_back( Order0Model(255) );
			
		_prevFieldPos.push_back(0);
		_currentFieldPos.push_back(0);
		_prevFieldValues.push_back(0);
		_currentFieldValues.push_back(0);
		_prevFieldTypes.push_back(FIELD_ASCII);
		_currentFieldTypes.push_back(FIELD_ASCII);
		_currentFieldZeroValues.push_back(0);
	}*/
	//cerr << "AbstractHeaderCoder::AbstractHeaderCoder - test seg 0.0" << endl;
	string _baseOutputname = _leon->_baseOutputname;
	string _baseInputname = _leon->_baseInputname;
	//cerr << "AbstractHeaderCoder::AbstractHeaderCoder - test seg 0.1" << endl;
	string typeModel[NB_MODELS];
	//cerr << "AbstractHeaderCoder::AbstractHeaderCoder - test seg 0.2" << endl;
	typeModel[TYPE_MODEL] = "_typeModel";
	typeModel[FIELD_INDEX_MODEL] = "_fieldIndexModel"; 
	typeModel[FIELD_COLUMN_MODEL] = "_fieldColumnModel"; 
	typeModel[MIS_SIZE_MODEL] = "_misSizeModel";
	typeModel[ASCII_MODEL] = "_asciiModel";
	typeModel[NUMERIC_HEAD_MODELS] = "_numericHeadModels";
	typeModel[ZERO_MODEL] = "_zeroModel";
	typeModel[HEADER_SIZE_MODEL] = "_headerSizeModel";
	//cerr << "AbstractHeaderCoder::AbstractHeaderCoder - test seg 0.3" << endl;
	//cerr << "AbstractHeaderCoder::AbstractHeaderCoder - test seg 1" << endl;

	if (_leon->_compress)
	{
		//cerr << "lol compression" << endl;
		//exit(EXIT_FAILURE);
		for(int i=0; i<NB_MODELS; ++i)
		{
			string ofstreamPath = _baseOutputname + ".noRangeEncoder" + typeModel[i];
			ofstreams[i].open(ofstreamPath, std::ofstream::app);
		}
	}
	if (_leon->_decompress)
	{
		//cerr << "lol _decompress" << endl;
		//exit(EXIT_FAILURE);
		for(int i=0; i<NB_MODELS; ++i)
		{
			string ifstreamPath = _baseInputname + ".noRangeEncoder" + typeModel[i];
			ifstreams[i].open(ifstreamPath, std::ifstream::in);
			//cerr << ifstreamPath << endl;
		}
		//exit(EXIT_FAILURE);
	}

	//cerr << "AbstractHeaderCoder::AbstractHeaderCoder - test seg 2" << endl;
}
	
void AbstractHeaderCoder::addFieldColumn(){
	
	_typeModel.push_back(Order0Model(HEADER_TYPE_COUNT+1));
	_fieldIndexModel.push_back(Order0Model(256));
	_fieldColumnModel.push_back(Order0Model(256));
	_misSizeModel.push_back(Order0Model(256));
	_asciiModel.push_back(Order0Model(128));
	_zeroModel.push_back(Order0Model(256));
	
	_numericModels.push_back(vector<Order0Model>());
	for(int j=0; j<CompressionUtils::NB_MODELS_PER_NUMERIC; j++)
		_numericModels[_numericModels.size()-1].push_back( Order0Model(256) );
		
	_prevFieldPos.push_back(0);
	_currentFieldPos.push_back(0);
	_prevFieldValues.push_back(0);
	_currentFieldValues.push_back(0);
	_prevFieldTypes.push_back(FIELD_ASCII);
	_currentFieldTypes.push_back(FIELD_ASCII);
	_prevFieldZeroValues.push_back(0);
	_currentFieldZeroValues.push_back(0);
}

int AbstractHeaderCoder::typeOfChar(u_int8_t c, bool* isDigit){
	if(isdigit(c)){
		*isDigit = true;
		return 1;
	}
	else if(isalpha(c)){
		*isDigit = false;
		return 1;
	}
	else{
		*isDigit = false;
		return 2;
	}
}

void AbstractHeaderCoder::splitHeader(){
	_fieldIndex = 0;
	_fieldStartPos = 0;
	_currentPos = 0;
	_isCurrentFieldNumeric = true;
	
	u_int8_t c;
	int charType;
	bool digitOnly;
	int lastCharType = typeOfChar(_currentHeader[0], &digitOnly);
	
	for(_currentPos=0; _currentPos<_currentHeader.size(); _currentPos++){
		c = _currentHeader[_currentPos];
		
		digitOnly = true;
		charType = typeOfChar(c, &digitOnly);
		
		if(charType != lastCharType){
			lastCharType = charType;
			makeField();
		}
		
		if(_isCurrentFieldNumeric){
			_isCurrentFieldNumeric = digitOnly;
		}
	}
	
	makeField();
	
	_currentFieldCount = _fieldIndex;
}

void AbstractHeaderCoder::makeField(){
	if(_fieldStartPos == _currentPos) return;
	
	//Adjust the maximum number fo field column
	int currentFieldColumn = _currentFieldPos.size();
	while(currentFieldColumn <= _fieldIndex+1){
		addFieldColumn();
		currentFieldColumn = _currentFieldPos.size();
	}
	
	_currentFieldPos[_fieldIndex] = _fieldStartPos;
	_currentFieldPos[_fieldIndex+1] = _currentPos;
	
	if(_isCurrentFieldNumeric){
		string field = _currentHeader.substr(_currentFieldPos[_fieldIndex], _currentFieldPos[_fieldIndex+1]-_currentFieldPos[_fieldIndex]);
		
		int zeroCount = 0;
		if(field[0] == '0'){
			while(field[0] == '0'){
				zeroCount += 1;
				field.erase(field.begin());
			}
		}
		
		//cout << " HAHAHAHAAHAHHAHA    " << field << endl; 
		_currentFieldZeroValues[_fieldIndex] = zeroCount;
		
		u_int64_t value = strtoul(field.c_str(), NULL, 0);
		_currentFieldValues[_fieldIndex] = value;
		
		if(zeroCount > 0){
			if(value == 0){
				_currentFieldTypes[_fieldIndex] = FIELD_ZERO_ONLY;
			}
			else{
				_currentFieldTypes[_fieldIndex] = FIELD_ZERO_AND_NUMERIC;
			}
		}
		else {
			_currentFieldTypes[_fieldIndex] = FIELD_NUMERIC;
		}
		//cout << "OOOOOOOOOOOOoo   " << strtoul("15313", NULL, 0) << endl;
	}
	else{
		_currentFieldTypes[_fieldIndex] = FIELD_ASCII;
	}
	
	#ifdef PRINT_DEBUG_ENCODER
		cout << "\tField "<< _fieldIndex << ": " << _currentHeader.substr(_currentFieldPos[_fieldIndex], _currentFieldPos[_fieldIndex+1]-_currentFieldPos[_fieldIndex]) << "  Digit only? " << _isCurrentFieldNumeric << endl;
	#endif
	_fieldIndex += 1;
	_fieldStartPos = _currentPos;
	_isCurrentFieldNumeric = true;
}

void AbstractHeaderCoder::endHeader(){
	_prevFieldCount = _currentFieldCount;
	
	#ifdef PRINT_DEBUG_ENCODER
		cout << "\tField count: " << _prevFieldCount << endl;
	#endif
	
	for(int i=0; i<_prevFieldCount+1; i++){
		_prevFieldPos[i] = _currentFieldPos[i];
		_prevFieldValues[i] = _currentFieldValues[i];
		_prevFieldTypes[i] = _currentFieldTypes[i];
		_prevFieldZeroValues[i] = _currentFieldZeroValues[i];
		
		_currentFieldZeroValues[i] = 0;
	}
	_prevHeader = _currentHeader;
	_misIndex = 0;
	_fieldIndex = 0;
	
	_processedSequenceCount += 1;
}

void AbstractHeaderCoder::startBlock(){

	_currentHeader = _leon->_firstHeader;
	cerr << "HeaderEncoder::operator - _currentHeader : " << _currentHeader<< endl;
	
	
	for(int i=0; i<_typeModel.size(); i++){
		_typeModel[i].clear();
		_fieldIndexModel[i].clear();
		_fieldColumnModel[i].clear();
		_misSizeModel[i].clear();
		_asciiModel[i].clear();
		_zeroModel[i].clear();
		
		for(int j=0; j<8; j++)
			_numericModels[i][j].clear();
			
	}
	_headerSizeModel.clear();
	
	splitHeader();
	endHeader();
	
	_processedSequenceCount = 0;
}

//====================================================================================
// ** HeaderEncoder
//====================================================================================
HeaderEncoder::HeaderEncoder(Leon* leon) :
AbstractHeaderCoder(leon) , _totalHeaderSize(0) ,_seqId(0)
{
	_thread_id = __sync_fetch_and_add (&_leon->_nb_thread_living, 1);

	
	//_firstHeader = firstHeader;
	//_rangeEncoder = new RangeEncoder();
}

HeaderEncoder::HeaderEncoder(const HeaderEncoder& copy) :
AbstractHeaderCoder(/*NULL*/copy._leon), _totalHeaderSize(0),_seqId(0)
{

	
	_leon = copy._leon;
	
	_thread_id = __sync_fetch_and_add (&_leon->_nb_thread_living, 1);
	startBlock();

	//_firstHeader = copy._firstHeader;
	//_rangeEncoder = new RangeEncoder();
}

HeaderEncoder::~HeaderEncoder(){
	
	
	if( _thread_id!=0 && (_seqId+1) % Leon::READ_PER_BLOCK != 0 ){
		writeBlock();
	}
	int nb_remaining = __sync_fetch_and_add (&_leon->_nb_thread_living, -1);
	
	__sync_fetch_and_add(&_leon->_totalHeaderSize, _totalHeaderSize);
	
#ifndef SERIAL
	//_leon->_blockwriter->incDone(1);
	_leon->_blockwriter->waitForWriter();
	
	if(nb_remaining==1)
	{
		_leon->_blockwriter->FlushWriter();
	}
#endif

	
	
	//if(_rangeEncoder->buffer.size() > 0){
		/*
		cout << "-----------------------" << _nb_occur<< endl;
		int thread_id = ((_lastSequenceIndex / Leon::READ_PER_BLOCK) % _leon->_nb_cores);
		cout << thread_id << endl;
		
		_rangeEncoder->flush();
		u_int8_t* buffer = &_rangeEncoder->buffer[0];
		_outputFile->fwrite(buffer, _rangeEncoder->buffer.size(), 1);
		
		cout << _lastSequenceIndex << endl;*/
	//}
	
	//delete _rangeEncoder;
}

/*
//static method
void HeaderEncoder::encodeFirstHeader(IFile* outputFile, const string& firstHeader){
	
	_rangeEncoder->encode(&_firstHeaderModel, firstHeader.size());
	for(int i=0; i<_currentHeader.size(); i++){
		_rangeEncoder->encode(&_firstHeaderModel, firstHeader[i]);
	}
	
	//splitHeader();
	//endHeader();
}*/

//int HeaderEncoder::getId(){
//	return ((_lastSequenceIndex / Leon::READ_PER_BLOCK) % _leon->_nb_cores);
//}

void HeaderEncoder::operator()(Sequence& sequence){

	_lastSequenceIndex = sequence.getIndex();
	_seqId = sequence.getIndex() ;


	_currentHeader = sequence.getComment();
	
	_totalHeaderSize += _currentHeader.size();


	processNextHeader();
	
	
	if(_processedSequenceCount >= Leon::READ_PER_BLOCK ){
		
		writeBlock();
		startBlock();
	}
	
}

void HeaderEncoder::writeBlock(){
	/*if(_rangeEncoder.getBufferSize() > 0){
		_rangeEncoder.flush();
	}*/
	
	int blockId = (  _seqId / Leon::READ_PER_BLOCK)   ;

	//printf("\nheader coder writeblock   bid %i   tid %i \n",blockId, _thread_id);
	cerr << "iofstreams_bufferSize : " << iofstreams_bufferSize << endl;
	//exit(EXIT_FAILURE);

	//_leon->writeBlock(_rangeEncoder.getBuffer(), _rangeEncoder.getBufferSize(), _processedSequenceCount, blockId);
	_leon->writeBlockNoRangeEncoder(iofstreams_bufferSize, _processedSequenceCount, blockId);
	//_rangeEncoder.clear();
}

void HeaderEncoder::processNextHeader(){
	
	#ifdef PRINT_DEBUG_ENCODER
		cout << _prevHeader << endl;
		cout << _currentHeader << endl;
	#endif
	splitHeader();
	compareHeader();
	endHeader();
}

void HeaderEncoder::compareHeader(){
	_fieldPos = 0;
	_misCurrentStartPos = -1;
	
	for(_fieldIndex=0; _fieldIndex<_currentFieldCount; _fieldIndex++){
		
		_currentFieldSize = _currentFieldPos[_fieldIndex+1] - _currentFieldPos[_fieldIndex];
		_prevFieldSize = _prevFieldPos[_fieldIndex+1]-_prevFieldPos[_fieldIndex];
		_misCurrentStartPos = -1;
		
		HeaderType prevFieldType = _prevFieldTypes[_fieldIndex];
		HeaderType currentFieldType = _currentFieldTypes[_fieldIndex];
		
		//Comparing numeric field
		if(prevFieldType == FIELD_NUMERIC && currentFieldType == FIELD_NUMERIC){
			#ifdef PRINT_DEBUG_ENCODER
				cout << "\t\tComparing numeric fields: " <<_prevFieldValues[_fieldIndex] << " " << _currentFieldValues[_fieldIndex] << endl;
			#endif
			if(_prevFieldValues[_fieldIndex] == _currentFieldValues[_fieldIndex]){ //match
				_lastMatchFieldIndex = _fieldIndex;
				continue;
			}
			//encodeNumeric();
		}
		//Comparing field with zero only
		else if(prevFieldType == FIELD_ZERO_ONLY && currentFieldType == FIELD_ZERO_ONLY){
			#ifdef PRINT_DEBUG_ENCODER
				cout << "\t\tComparing fields with zero only: "  << endl;
			#endif
			if(_prevFieldZeroValues[_fieldIndex] == _currentFieldZeroValues[_fieldIndex]){ //match
				_lastMatchFieldIndex = _fieldIndex;
				continue;
			}
			//encodeNumeric();
		}
		//Comparing field with zero at begining and numeric
		else if(prevFieldType == FIELD_ZERO_AND_NUMERIC && currentFieldType == FIELD_ZERO_AND_NUMERIC){
			#ifdef PRINT_DEBUG_ENCODER
				cout << "\t\tComparing fields with zero at begining and numeric: "  << endl;
			#endif
			if(_prevFieldZeroValues[_fieldIndex] == _currentFieldZeroValues[_fieldIndex] && _prevFieldValues[_fieldIndex] == _currentFieldValues[_fieldIndex]){ //match
				_lastMatchFieldIndex = _fieldIndex;
				continue;
			}
			//encodeNumeric();
		}
		
		
		//Encoding numeric field
		if(currentFieldType == FIELD_NUMERIC || currentFieldType == FIELD_ZERO_ONLY || currentFieldType == FIELD_ZERO_AND_NUMERIC){
			encodeNumeric();
		}
		//Comparing ascii fields
		else{
			for(_fieldPos=0; _fieldPos<_currentFieldSize; _fieldPos++){
			
				if(_fieldIndex >= _prevFieldCount){
					_misCurrentStartPos = _fieldPos;
					break;
				}
					
				if(_fieldPos >= _prevFieldSize){
					_misCurrentStartPos = _fieldPos;
					break;
				}

				u_int8_t c = _currentHeader[_currentFieldPos[_fieldIndex]+_fieldPos];
				
				#ifdef PRINT_DEBUG_ENCODER
					cout << "\t\tComparing: " << _prevHeader[_prevFieldPos[_fieldIndex]+_fieldPos] << " " << c << "  ";
				#endif
				if(c == _prevHeader[_prevFieldPos[_fieldIndex]+_fieldPos]){
					#ifdef PRINT_DEBUG_ENCODER
						cout << "match" << endl;
					#endif
				}
				else{
					#ifdef PRINT_DEBUG_ENCODER
						cout << "mismatch" << endl;
					#endif
					_misCurrentStartPos = _fieldPos;
					break;
				}
			}
			
			if(_misCurrentStartPos != -1){
				encodeAscii();
			}
			else if(_fieldPos != _prevFieldSize){ //All the character of the current field match but there are always character in the current field of the prev header
				_misCurrentStartPos = _fieldPos;
				encodeAscii();
			}
			else{
				_lastMatchFieldIndex = _fieldIndex;
			}
			
		}
		
	}
	
	//if(_currentFieldPos[_fieldIndex]+_fieldPos == _currentHeader.size()){
	//	_misCurrentStartPos = _fieldPos;
	//	encodeMismatch();
	//}
	
	//cout << _lastMatchFieldIndex << " " << _fieldIndex << endl;
	
	//if the last field match, we have to signal to the decoder to add the last matching field of the prev header
	
	if(_lastMatchFieldIndex == _fieldIndex-1){
		//_rangeEncoder.encode(_typeModel[_misIndex], HEADER_END_MATCH);
		ofstreams[TYPE_MODEL] << (int) HEADER_END_MATCH << endl;
		//cerr << "HEADER_END_MATCH : " << HEADER_END_MATCH << endl;
		//cerr << "test : " << std::to_string((int)HEADER_END_MATCH).length() << endl;
		//exit(EXIT_FAILURE);
		iofstreams_bufferSize += std::to_string((int) HEADER_END_MATCH).length();
		//_rangeEncoder.encode(_headerSizeModel, _currentHeader.size());
		ofstreams[HEADER_SIZE_MODEL] << (int) _currentHeader.size() << endl;
		iofstreams_bufferSize += std::to_string((int) _currentHeader.size()).length();
		//_misCurrentStartPos = _currentFieldSize;
		//encodeAscii();
	}
	else{
		//_rangeEncoder.encode(_typeModel[_misIndex], HEADER_END);
		ofstreams[TYPE_MODEL] << (int) HEADER_END << endl;
		iofstreams_bufferSize += std::to_string((int) HEADER_END).length();
	}
	//_misIndex += 1;
	
	//end of header
	//_rangeEncoder.encode(_typeModel[_misIndex], HEADER_END);
	//_rangeEncoder.encode(_headerSizeModel, _currentHeader.size());
	
}





/*
void HeaderEncoder::encode(){
	if(_misCurrentStartPos != -1){
		encodeMismatch();
	}
	else if(_fieldPos != _prevFieldSize){ //All the character of the current field match but there are always character in the current field of the prev header
		_misCurrentStartPos = _fieldPos;
		encodeMismatch();
	}
	else if(_currentFieldPos[_fieldIndex]+_fieldPos == _currentHeader.size()){
		_misCurrentStartPos = _fieldPos;
		encodeMismatch();
	}
}
void HeaderEncoder::encodeMismatch(){
	HeaderType fieldType = _currentFieldTypes[_fieldIndex];
	if(fieldType == FIELD_NUMERIC || fieldType == FIELD_ZERO_ONLY || fieldType == FIELD_ZERO_AND_NUMERIC){
		encodeNumeric();
	}
	else{
		encodeAscii();
	}
	_misIndex += 1;
}*/

void HeaderEncoder::encodeNumeric(){
	u_int64_t zeroCount = _currentFieldZeroValues[_fieldIndex];
	u_int64_t fieldValue = _currentFieldValues[_fieldIndex];
	
	HeaderType currentFieldType = _currentFieldTypes[_fieldIndex];
	
	if(currentFieldType == FIELD_ZERO_ONLY){
		#ifdef PRINT_DEBUG_ENCODER
			cout << "\t\t\tField with zero only" << endl;
			cout << "\t\t\tEnconding zero count: " << zeroCount << endl;
		#endif
		//_rangeEncoder.encode(_typeModel[_misIndex], FIELD_ZERO_ONLY);
		ofstreams[TYPE_MODEL] << (int) FIELD_ZERO_ONLY << endl;
		iofstreams_bufferSize += std::to_string((int) FIELD_ZERO_ONLY).length();
		//_rangeEncoder.encode(_fieldIndexModel[_misIndex], _fieldIndex);
		ofstreams[FIELD_INDEX_MODEL] << (int) _fieldIndex << endl;
		iofstreams_bufferSize += std::to_string((int) _fieldIndex).length();
		//_rangeEncoder.encode(_zeroModel[_misIndex], zeroCount);
		ofstreams[ZERO_MODEL] << (int) zeroCount << endl;
		iofstreams_bufferSize += std::to_string((int) zeroCount).length();
		_misIndex += 1;
		return;
	}
	else if(currentFieldType == FIELD_ZERO_AND_NUMERIC){
		#ifdef PRINT_DEBUG_ENCODER
			cout << "\t\t\tField with zero and numeric" << endl;
			cout << "\t\t\tEnconding zero count: " << zeroCount << endl;
		#endif
		//_rangeEncoder.encode(_typeModel[_misIndex], FIELD_ZERO_AND_NUMERIC);
		ofstreams[TYPE_MODEL] << (int) FIELD_ZERO_AND_NUMERIC << endl;
		iofstreams_bufferSize += std::to_string((int) FIELD_ZERO_AND_NUMERIC).length();
		//_rangeEncoder.encode(_fieldIndexModel[_misIndex], _fieldIndex);
		ofstreams[FIELD_INDEX_MODEL] << (int) _fieldIndex << endl;
		iofstreams_bufferSize += std::to_string((int) _fieldIndex).length();
		//_rangeEncoder.encode(_zeroModel[_misIndex], zeroCount);
		ofstreams[ZERO_MODEL] << (int) zeroCount << endl;
		iofstreams_bufferSize += std::to_string((int) zeroCount).length();
		_misIndex += 1;
	}
	
	/*
	if(zeroCount > 0){
		if(fieldValue == 0){
			#ifdef PRINT_DEBUG_ENCODER
				cout << "\t\t\tField with zero only" << endl;
			#endif
			_rangeEncoder.encode(_typeModel[_misIndex], FIELD_ZERO_ONLY);
		}
		else{
			#ifdef PRINT_DEBUG_ENCODER
				cout << "\t\t\tField with zero at begining and numeric" << endl;
			#endif
			_rangeEncoder.encode(_typeModel[_misIndex], FIELD_ZERO_AND_NUMERIC);
		}
		_rangeEncoder.encode(_fieldIndexModel[_misIndex], _fieldIndex);
		//cout << "sdf " << _fieldIndex << endl;
		_rangeEncoder.encode(_zeroModel[_misIndex], zeroCount);
		//#ifdef PRINT_DEBUG_ENCODER
		//	cout << "\t\t\t\tZero count: " << zeroCount << endl;
		//#endif
		//_currentFieldZeroValues[_fieldIndex] = 0;
		//_misIndex += 1;
		//if(fieldValue == 0) return;
	}*/
	
	//cout << _currentFieldPos[_fieldIndex] << endl;
	//cout << _currentPos << endl;
	//u_int64_t fieldValue = strtoul(_currentHeader.substr(_currentFieldPos[_fieldIndex], _currentFieldPos[_fieldIndex]-_currentPos).c_str(), NULL, 0);
	
	//if(fieldValue == 0){
	//	_misIndex -= 1;
	//	return;
	//}
	u_int64_t value = fieldValue;
	u_int64_t prevValue = _prevFieldValues[_fieldIndex];


	//int valueByteCount = CompressionUtils::getByteCount(value);
	
	#ifdef PRINT_DEBUG_ENCODER
		cout << "\t\t\tPrev value: " << prevValue << endl;
		cout << "\t\t\tField value: " << value << "    Byte: " << valueByteCount << endl;
	#endif
	
	u_int64_t deltaValue;
	int deltaType = CompressionUtils::getDeltaValue(value, prevValue, &deltaValue);
	
	if(deltaType == 0){
		//_rangeEncoder.encode(_typeModel[_misIndex], FIELD_NUMERIC);
		ofstreams[TYPE_MODEL] << (int) FIELD_NUMERIC << endl;
		iofstreams_bufferSize += std::to_string((int) FIELD_NUMERIC).length();
	}
	else if(deltaType == 1){
		//_rangeEncoder.encode(_typeModel[_misIndex], FIELD_DELTA);
		ofstreams[TYPE_MODEL] << (int) FIELD_DELTA << endl;
		iofstreams_bufferSize += std::to_string((int) FIELD_DELTA).length();
		value = deltaValue;
	}
	else if(deltaType == 2){
		//_rangeEncoder.encode(_typeModel[_misIndex], FIELD_DELTA_2);
		ofstreams[TYPE_MODEL] << (int) FIELD_DELTA_2 << endl;
		iofstreams_bufferSize += std::to_string((int) FIELD_DELTA_2).length();
		value = deltaValue;
	}
	/*
	//if(prevValue <= value){
	u_int64_t deltaValue = value - prevValue;
	int deltaByteCount = CompressionUtils::getByteCount(deltaValue);
	
	#ifdef PRINT_DEBUG_ENCODER
		cout << "\t\t\tDelta value: " << deltaValue << "    Byte: " << deltaByteCount << endl;
	#endif
	if(deltaValue >= 0 && deltaByteCount <= valueByteCount){
		_rangeEncoder.encode(_typeModel[_misIndex], FIELD_DELTA);
		value = deltaValue;
		//valueByteCount = deltaByteCount;
	}
	else{
	
		
		u_int64_t deltaValue2 = prevValue - value;
		int deltaByteCount2 = CompressionUtils::getByteCount(deltaValue2);
	
		if(deltaValue2 >= 0 && deltaByteCount2 <= valueByteCount){
			_rangeEncoder.encode(_typeModel[_misIndex], FIELD_DELTA_2);
			value = deltaValue2;
			//valueByteCount = deltaByteCount2;
		}
		else{
			_rangeEncoder.encode(_typeModel[_misIndex], FIELD_NUMERIC);
		}
	}
	//}
	//else{
	//	_rangeEncoder.encode(_typeModel[_misIndex], FIELD_NUMERIC);
	//}
	*/

		
	  
	//_rangeEncoder.encode(_fieldIndexModel[_misIndex], _fieldIndex);
	ofstreams[FIELD_INDEX_MODEL] << (int) _fieldIndex << endl;
	iofstreams_bufferSize += std::to_string((int) _fieldIndex).length();
	//CompressionUtils::encodeNumeric(_rangeEncoder, _numericModels[_misIndex], value);
	//cerr << "value : " << value<< endl;
	//exit(EXIT_FAILURE);
	ofstreams[NUMERIC_HEAD_MODELS] << (int) value << endl;
	iofstreams_bufferSize += std::to_string((int) value).length();
	//_rangeEncoder->encode(&_fieldColumnModel[_misIndex], 0);
	//_prevFieldValues[_fieldIndex] = fieldValue;
	
	_misIndex += 1;
}

void HeaderEncoder::encodeAscii(){
	int missSize = _currentFieldSize - _misCurrentStartPos;//_currentPos - _misCurrentStartPos;
	//cout << _currentFieldSize <<  " " << _fieldPos << endl;
	//_rangeEncoder.encode(_typeModel[_misIndex], FIELD_ASCII);
	ofstreams[TYPE_MODEL] << (int) FIELD_ASCII << endl;
	iofstreams_bufferSize += std::to_string((int) FIELD_ASCII).length();
	//_rangeEncoder.encode(_fieldIndexModel[_misIndex], _fieldIndex);
	ofstreams[FIELD_INDEX_MODEL] << (int) _fieldIndex << endl;
	iofstreams_bufferSize += std::to_string((int) _fieldIndex).length();
	//_rangeEncoder.encode(_fieldColumnModel[_misIndex], _misCurrentStartPos);
	ofstreams[FIELD_COLUMN_MODEL] << (int) _misCurrentStartPos << endl;
	iofstreams_bufferSize += std::to_string((int) _misCurrentStartPos).length();
	//_rangeEncoder.encode(_misSizeModel[_misIndex], missSize);
	ofstreams[MIS_SIZE_MODEL] << (int) missSize << endl;
	iofstreams_bufferSize += std::to_string((int) missSize).length();
	#ifdef PRINT_DEBUG_ENCODER
		cout << "\t\t\t<Mismatch> " << "    Type: " << "ASCII" << "    Field: " << _fieldIndex << "    Column: " << _misCurrentStartPos << "    Size: " << missSize << endl;
	#endif
	//for(int j=_misCurrentStartPos; j<_currentPos; j++){
	for(int i=_misCurrentStartPos; i < _misCurrentStartPos+missSize; i++){
		#ifdef PRINT_DEBUG_ENCODER
			cout << "\t\t\tEncoding: " << _currentHeader[_currentFieldPos[_fieldIndex]+i] << endl;
		#endif
		//cout << _currentHeader[j] << flush;
		//_rangeEncoder.encode(_asciiModel[_misIndex], _currentHeader[_currentFieldPos[_fieldIndex]+i]);
		ofstreams[ASCII_MODEL] << (int) _currentHeader[_currentFieldPos[_fieldIndex]+i] << endl;
		iofstreams_bufferSize += std::to_string((int) _currentHeader[_currentFieldPos[_fieldIndex]+i]).length();
	}

	_misIndex += 1;
}











//====================================================================================
// ** HeaderDecoder
//====================================================================================
HeaderDecoder::HeaderDecoder(Leon* leon, const string& inputFilename) :
AbstractHeaderCoder(leon)
//, _rangeDecoder(inputFile)
{
	_inputFile = new ifstream(inputFilename.c_str(), ios::in|ios::binary);
	_finished = false;
	//_outputFile = outputFile;
	// = new RangeDecoder(inputFile);
	
	//_outputFile = new OutputFile(outputFile);

	
}

HeaderDecoder::HeaderDecoder(Leon* leon, Requests* req, const string& inputFilename) :
AbstractHeaderCoder(leon)
{
	//cout << endl << "debug - HeaderDecoder - constructor2 - inputFilename : " << inputFilename << endl;
	_rangeDecoder = req->_rangeDecoder;
	_inputFile = new ifstream(inputFilename.c_str(), ios::in|ios::binary);
	_finished = false;
	//_outputFile = outputFile;
	// = new RangeDecoder(inputFile);
	
	//_outputFile = new OutputFile(outputFile);
	
}

HeaderDecoder::~HeaderDecoder(){
	//delete _rangeDecoder;
	//delete _outputFile;
	delete _inputFile;
}

void HeaderDecoder::setup(u_int64_t blockStartPos, u_int64_t blockSize, int sequenceCount){
	
	//cout << endl << "debug - HeaderDecoder - setup - blockStartPos : " << blockStartPos << endl;
	//cout << "debug - HeaderDecoder - setup - blockSize : " << blockSize << endl;
	//cout << "debug - HeaderDecoder - setup - sequenceCount : " << sequenceCount << endl;



	startBlock();
	_rangeDecoder.clear();
	
	_inputFile->seekg(blockStartPos, _inputFile->beg);
	_rangeDecoder.setInputFile(_inputFile);
	
	_blockStartPos = blockStartPos;
	_blockSize = blockSize;
	
	#ifdef PRINT_DEBUG_DECODER
		cout << "\t-----------------------" << endl;
		cout << "\tDecoding block " << _blockStartPos << " - " << _blockStartPos+_blockSize << endl;
	#else
		//_leon->_progress_decode->inc(1);

	#endif
	
	//_currentHeader = _leon->_firstHeader;
	//endHeader();
	_currentHeader.clear();
	_misIndex = 0;
	
	_sequenceCount = sequenceCount;
	
}

void HeaderDecoder::execute(){
	//cerr << endl << "DEBUG HEADER CODER EXECUTE BEGIN" << endl;
	//cout << "executing" << endl;
	//decodeFirstHeader();
	cerr << "HeaderDecoder::execute - BEGIN" << endl;
	cerr << "HeaderDecoder::execute - _processedSequenceCount : " << _processedSequenceCount << endl;
	cerr << "HeaderDecoder::execute - _sequenceCount : " << _sequenceCount << endl;
	while(_processedSequenceCount < _sequenceCount){
		//cerr << "HeaderDecoder::execute - HERE" << endl;
		//int i=0;
		//while(_inputFile->tellg() <= _blockStartPos+_blockSize){
			
			//cout << "tellg: " << _inputFile->tellg() << endl;
			//if(i>= 1) return;
		//for(u_int64_t i=0; i<_blockSize; i++){
			
		//while(!_inputFile->eof()){
		//cerr << "debug header decoder : readType before" << endl;
		int type;// = _rangeDecoder.nextByte(_typeModel[_misIndex]);
		//std::getline(noRangeDecoderIfstream_numericModel, _textDecodeInfo);
		//std::getline(ifstreams[TYPE_MODEL], _textDecodeInfo);
		//std::getline(ifstreams[TYPE_MODEL], _textDecodeInfo);
		//cerr << _textDecodeInfo << endl;
		//int type2;
		ifstreams[TYPE_MODEL] >> type;
		//ifstreams[TYPE_MODEL] >> type2;
		cerr << "isopen : " << ifstreams[TYPE_MODEL].is_open() << endl;
		cerr << "TYPE_MODEL : " << type << endl;
				//exit(EXIT_FAILURE);
		//cerr << "debug header decoder : readType after" << endl;
		//cerr << "debug header decoder : type = " << (int) type << endl;
		#ifdef PRINT_DEBUG_DECODER
			cout << "\t\tNext type is: " << (int) type << endl;
		#endif
		//if(_inputFile->eof()) return;
		//if(type==0){
			//return; //Pourquoi y a t-il des zero a la fin du fichier ????????
		//}
		//cout << "type again " << (int)type << endl;

		
		if(type == HEADER_END){
			/*
			u_int64_t headerSize = CompressionUtils::decodeNumeric(_rangeDecoder, _numericSizeModel[_misIndex], _numericModels[_misIndex]);
			int i = _prevFieldPos[_fieldIndex];
			while(_currentHeader.size() < headerSize){
				//_currentHeader += _prevHeader
			}*/
	//		cout << "debug header decoder : enter endHeader" << endl;
			endHeader();
	//		cout << "debug header decoder : exit endHeader" << endl;
			//i+=1;
		}
		else if(type == HEADER_END_MATCH){
			//decodeMatch();
			//cout << "debug header decoder : headerSize before" << endl;
			/*u_int8_t*/int headerSize;// = _rangeDecoder.nextByte(_headerSizeModel);
			//std::getline(ifstreams[HEADER_SIZE_MODEL], _textDecodeInfo);
			//std::getline(ifstreams[HEADER_SIZE_MODEL], _textDecodeInfo);
			//cerr << _textDecodeInfo << endl;
			ifstreams[HEADER_SIZE_MODEL] >> headerSize;
			cerr << "HEADER_SIZE_MODEL : " << headerSize << endl;
			//cout << "debug header decoder : headerSize before" << endl;

			for(/*_fieldIndex*/; _fieldIndex < _prevFieldCount; _fieldIndex++){
				#ifdef PRINT_DEBUG_DECODER
					cout << "\t\t\tAdding from prev header: " << _prevHeader.substr(_prevFieldPos[_fieldIndex], _prevFieldPos[_fieldIndex+1]-_prevFieldPos[_fieldIndex]) << endl;
				#endif
				_currentHeader += _prevHeader.substr(_prevFieldPos[_fieldIndex], _prevFieldPos[_fieldIndex+1]-_prevFieldPos[_fieldIndex]);
				if(_currentHeader.size() >= headerSize) break;
			}
			cerr << "_currentHeader : " << _currentHeader << endl;	
			
			endHeader();
		}
		else{
	//		cout << "debug header decoder : enter decodeMatch" << endl;
			decodeMatch();
	//		cout << "debug header decoder : enter decodeMatch" << endl;

			if(type == FIELD_ASCII){
				decodeAscii();
				_fieldIndex += 1;
				_misIndex += 1;
			}
			else if(type == FIELD_NUMERIC){
				decodeNumeric();
				_fieldIndex += 1;
				_misIndex += 1;
			}
			else if(type == FIELD_DELTA){
				decodeDelta();
				_fieldIndex += 1;
				_misIndex += 1;
			}
			else if(type == FIELD_DELTA_2){
				decodeDelta2();
				_fieldIndex += 1;
				_misIndex += 1;
			}
			else if(type == FIELD_ZERO_ONLY){
				decodeZero();
				_fieldIndex += 1;
				_misIndex += 1;
			}
			else if(type == FIELD_ZERO_AND_NUMERIC){
				decodeZero();
				_misIndex += 1;
				//decodeNumeric();
				//_fieldIndex += 1;
			}
			//_prevPos = _prevFieldPos[_fieldIndex+1];

		}
		
		//cout << "lala" << endl;
	}
	
	_finished = true;
}
/*
void HeaderDecoder::decodeFirstHeader(){
	u_int8_t size = _rangeDecoder->nextByte(&_firstHeaderModel);
	for(int i=0; i < size; i++){
		u_int8_t c = _rangeDecoder->nextByte(&_firstHeaderModel);
		_currentHeader += c;
	}
	endHeader();
}
*/

void HeaderDecoder::decodeMatch(){
	//cout << "debug header decoder : decodeMAtch read misFieldIndex before" << endl;
	/*u_int8_t*/ int misFieldIndex;// = _rangeDecoder.nextByte(_fieldIndexModel[_misIndex]);
	//std::getline(ifstreams[FIELD_INDEX_MODEL], _textDecodeInfo);
	// std::getline(ifstreams[FIELD_INDEX_MODEL], _textDecodeInfo);
	// cerr << _textDecodeInfo << endl;
	ifstreams[FIELD_INDEX_MODEL] >> misFieldIndex;
	cerr << "FIELD_INDEX_MODEL : " << misFieldIndex << endl;
	//cout << "debug header decoder : decodeMAtch read misFieldIndex after" << endl;
	//cout << "debug header decoder : misFieldIndex = " << (int) misFieldIndex << endl;
	#ifdef PRINT_DEBUG_DECODER
		cout << "\t\tMatch to field: " << (int)misFieldIndex << endl;
	#endif
	//				cerr << "HERE" << endl;
	cerr << "_fieldIndex : " << (int) _fieldIndex << endl;
	cerr << "_prevHeader : " << _prevHeader << endl;
	for(/*_fieldIndex*/; _fieldIndex < misFieldIndex; _fieldIndex++){
		//#ifdef PRINT_DEBUG_DECODER
		//	cout << "\t\t\tAdding from prev header: " << _prevHeader.substr(_prevFieldPos[_fieldIndex], _prevFieldPos[_fieldIndex+1]-_prevFieldPos[_fieldIndex]) << endl;
		//#endif
		cerr << "_fieldIndex : " << (int) _fieldIndex << endl;
		cerr << "misFieldIndex : " << (int) misFieldIndex << endl;
		_currentHeader += _prevHeader.substr(_prevFieldPos[_fieldIndex], _prevFieldPos[_fieldIndex+1]-_prevFieldPos[_fieldIndex]);
		//cerr << "_currentHeader : " << _currentHeader << endl;
	}
	cerr << "_currentHeader : " << _currentHeader << endl;	
	//cerr << "HERE" << endl;
			//exit(EXIT_FAILURE);
}

void HeaderDecoder::decodeAscii(){
	/*u_int8_t*/int misColumn;// = _rangeDecoder.nextByte(_fieldColumnModel[_misIndex]);
	// std::getline(ifstreams[FIELD_COLUMN_MODEL], _textDecodeInfo);
	// std::getline(ifstreams[FIELD_COLUMN_MODEL], _textDecodeInfo);
	// cerr << _textDecodeInfo << endl;
	ifstreams[FIELD_COLUMN_MODEL] >> misColumn;
	cerr << "FIELD_COLUMN_MODEL : " << misColumn << endl;

	/*u_int8_t*/int misSize;// = _rangeDecoder.nextByte(_misSizeModel[_misIndex]);
	// std::getline(ifstreams[MIS_SIZE_MODEL], _textDecodeInfo);
	// std::getline(ifstreams[MIS_SIZE_MODEL], _textDecodeInfo);
	// cerr << _textDecodeInfo << endl;
	ifstreams[MIS_SIZE_MODEL] >> misSize;
	cerr << "MIS_SIZE_MODEL : " << misSize << endl;
	#ifdef PRINT_DEBUG_DECODER
		cout << "\t\tDecoding   Type: ASCII     Column: " << (int)misColumn << "    Size: " << (int)misSize << endl;
	#endif
	
	if(_fieldIndex < _prevFieldCount){
		for(int fieldPos=0; fieldPos<misColumn; fieldPos++){
			_currentHeader += _prevHeader[_prevFieldPos[_fieldIndex]+fieldPos];
		}
		cerr << "_currentHeader : " << _currentHeader << endl;	
	}
	
	for(int i=0; i<misSize; i++){
		/*u_int8_t*/int c;// = _rangeDecoder.nextByte(_asciiModel[_misIndex]);
		// std::getline(ifstreams[ASCII_MODEL], _textDecodeInfo);
		// std::getline(ifstreams[ASCII_MODEL], _textDecodeInfo);
		// cerr << _textDecodeInfo << endl;
		ifstreams[ASCII_MODEL] >> c;
		cerr << "ASCII_MODEL : " << c << endl;
		#ifdef PRINT_DEBUG_DECODER
			cout << "\t\t\tAdding: " << c << " (" << (int)c << ")"<< endl;
		#endif
		//_currentHeader2[_currentPos] = c;
		_currentHeader += c;
		//_currentPos += 1;
	}
	cerr << "_currentHeader : " << _currentHeader << endl;	
	
}	

void HeaderDecoder::decodeNumeric(){
	//u_int8_t misSize = _rangeDecoder.nextByte(_misSizeModel[_misIndex]);
	#ifdef PRINT_DEBUG_DECODER
		cout << "\t\tDecoding   Type: NUMERIC" << endl; //"    Size: " << (int)misSize << endl;
	#endif
	
	int value;// = CompressionUtils::decodeNumeric(_rangeDecoder, _numericModels[_misIndex]);
	// std::getline(ifstreams[NUMERIC_HEAD_MODELS], _textDecodeInfo);
	// std::getline(ifstreams[NUMERIC_HEAD_MODELS], _textDecodeInfo);
	// cerr << _textDecodeInfo << endl;
	ifstreams[NUMERIC_HEAD_MODELS] >> value;
	cerr << "NUMERIC_HEAD_MODELS : " << value << endl;
	//_currentHeader += CompressionUtils::numberToString(value);
	
	char temp[200];
	snprintf(temp,200,"%llu",value);
	_currentHeader += string(temp);
	//_currentHeader += to_string(value); // C++11
	cerr << "_currentHeader : " << _currentHeader << endl;	

	#ifdef PRINT_DEBUG_DECODER
		cout << "\t\t\tAdding: " << string(temp) << endl;
	#endif
}

void HeaderDecoder::decodeDelta(){
	//u_int8_t misSize = _rangeDecoder.nextByte(_misSizeModel[_misIndex]);
	#ifdef PRINT_DEBUG_DECODER
		cout << "\t\tDecoding   Type: DELTA" << endl;//"    Size: " << (int)misSize << endl;
	#endif
	
	u_int64_t value;// = CompressionUtils::decodeNumeric(_rangeDecoder, _numericModels[_misIndex]);
	// std::getline(ifstreams[NUMERIC_HEAD_MODELS], _textDecodeInfo);
	// std::getline(ifstreams[NUMERIC_HEAD_MODELS], _textDecodeInfo);
	// cerr << _textDecodeInfo << endl;
	ifstreams[NUMERIC_HEAD_MODELS] >> value;
	cerr << "NUMERIC_HEAD_MODELS : " << value << endl;
	/*
	u_int64_t value = 0;
	for(int i=0; i<misSize; i++){
		u_int8_t byteValue = _rangeDecoder.nextByte(_numericModels[_misIndex][i]);
		value |= (byteValue << i*8);
	}
	cout << "lala  " << _prevFieldValues[_fieldIndex] << endl;*/
	//value += _prevFieldValues[_fieldIndex];
	value = CompressionUtils::getValueFromDelta(1, _prevFieldValues[_fieldIndex], value);
	
	char temp[200];
	snprintf(temp,200,"%llu",value);
	_currentHeader += string(temp);
	//_currentHeader += to_string(value);
	cerr << "_currentHeader : " << _currentHeader << endl;	
	#ifdef PRINT_DEBUG_DECODER
		cout << "\t\t\tAdding: " << string(temp) << endl;
	#endif
}

void HeaderDecoder::decodeDelta2(){
	//u_int8_t misSize = _rangeDecoder.nextByte(_misSizeModel[_misIndex]);
	#ifdef PRINT_DEBUG_DECODER
		cout << "\t\tDecoding   Type: DELTA 2" << endl;//"    Size: " << (int)misSize << endl;
	#endif
	
	u_int64_t value;// = CompressionUtils::decodeNumeric(_rangeDecoder, _numericModels[_misIndex]);
	// std::getline(ifstreams[NUMERIC_HEAD_MODELS], _textDecodeInfo);
	// std::getline(ifstreams[NUMERIC_HEAD_MODELS], _textDecodeInfo);
	// cerr << _textDecodeInfo << endl;
	ifstreams[NUMERIC_HEAD_MODELS] >> value;
	cerr << "NUMERIC_HEAD_MODELS : " << value << endl;
	/*
	u_int64_t value = 0;
	for(int i=0; i<misSize; i++){
		u_int8_t byteValue = _rangeDecoder.nextByte(_numericModels[_misIndex][i]);
		value |= (byteValue << i*8);
	}
	cout << "lala  " << _prevFieldValues[_fieldIndex] << endl;*/
	//value = _prevFieldValues[_fieldIndex] - value;
	value = CompressionUtils::getValueFromDelta(2, _prevFieldValues[_fieldIndex], value);
	char temp[200];
	snprintf(temp,200,"%llu",value);
	_currentHeader += string(temp);
	//_currentHeader += to_string(value);
	cerr << "_currentHeader : " << _currentHeader << endl;	
	#ifdef PRINT_DEBUG_DECODER
		cout << "\t\t\tAdding: " << string(temp) << endl;
	#endif
}

void HeaderDecoder::decodeZero(){
	/*u_int8_t*/int zeroCount;// = _rangeDecoder.nextByte(_zeroModel[_misIndex]);
	// std::getline(ifstreams[ZERO_MODEL], _textDecodeInfo);
	// std::getline(ifstreams[ZERO_MODEL], _textDecodeInfo);
	// cerr << _textDecodeInfo << endl;
	ifstreams[ZERO_MODEL] >> zeroCount;
	cerr << "ZERO_MODEL : " << zeroCount << endl;

	#ifdef PRINT_DEBUG_DECODER
		cout << "\t\tDecoding   Type: ZERO     Size: " << (int)zeroCount << endl;
	#endif
	
	for(int i=0; i<zeroCount; i++){
		
		#ifdef PRINT_DEBUG_DECODER
			cout << "\t\t\tAdding: 0"<< endl;
		#endif
		
		_currentHeader += '0';
		cerr << "_currentHeader : " << _currentHeader << endl;	
	}
}

void HeaderDecoder::endHeader(){
	_buffer += _currentHeader + '\n';
	//_outputFile->write((_currentHeader+'\n').c_str(), _currentHeader.size()+1);
	
	//_currentHeader.erase(_currentHeader.begin());
	//cout << _currentHeader << endl;
	//for(int i=0; i<_currentHeader.size(); i++){
		
		//_outputFile->writeByte(_currentHeader[i]);
		//_prevHeader2[i] = _currentHeader[i];
	//}
	//_outputFile->writeByte('\n');
	
	#ifdef PRINT_DEBUG_DECODER
		cout << _currentHeader << endl;
		//for(int i=0; i<_currentPos; i++){
		//	cout << _currentHeader2[i];
		//}
	#endif
	
	
	splitHeader();
	
	
	//_prevHeaderSize = _currentPos;
	//_prevFieldCount = _fieldIndex;
	
	AbstractHeaderCoder::endHeader();
	//for(int i=0; i<_prevFieldCount+1; i++){
	//	_prevFieldPos[i] = _currentFieldPos[i];
	//}
		
		
	//_prevHeader = _currentHeader;
	_currentHeader.clear();
	_misIndex = 0;
	//_prevPos = 0;
}


