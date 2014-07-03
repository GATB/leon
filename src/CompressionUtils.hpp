
#ifndef _COMPRESSIONUTILS_HPP_
#define _COMPRESSIONUTILS_HPP_

using namespace std;
class Leon;
//====================================================================================
// ** COMPRESSIONUTILS
//====================================================================================
class CompressionUtils
{
	public:
	
		
		
		//Take a value and return its number of byte
		static int getByteCount(u_int64_t n){
			if(n < 0x100){
				return 1;
			}
			else if(n >= 0x100 && n < 0x10000){
				return 2;
			}
			else if(n >= 0x10000 && n < 0x1000000){
				return 3;
			}
			else if(n >= 0x1000000 && n < 0x100000000){
				return 4;
			}
			else if(n >= 0x100000000 && n < 0x10000000000){
				return 5;
			}
			else if(n >= 0x10000000000 && n < 0x1000000000000){
				return 6;
			}
			else if(n >= 0x1000000000000 && n < 0x100000000000000){
				return 7;
			}
			else{
				return 8;
			}

		}
		
		//Encode a numeric value
		//First encode the byte count of the value
		//Then encode the value of each byte
		static void encodeNumeric(RangeEncoder& rangeEncoder, Order0Model& byteCountModel, vector<Order0Model>& numericModels, u_int64_t value){
			int valueByteCount = getByteCount(value);
			//cout << "Utils: " << (int)valueByteCount << endl;
			rangeEncoder.encode(byteCountModel, valueByteCount);
				
			for(int i=0; i<valueByteCount; i++){
				//cout << "Utils: " << ((value >> i*8) & 0xff) << endl;
				rangeEncoder.encode(numericModels[i], (value >> i*8) & 0xff);
			}
		}

		static void encodeFixedNumeric(RangeEncoder& rangeEncoder, vector<Order0Model>& numericModels, u_int64_t value, int byteCount){
			//cout << "sdf" << value << " " << byteCount << endl;
			for(int i=0; i<byteCount; i++){
				//cout << ((value >> i*8) & 0xff) << endl;
				rangeEncoder.encode(numericModels[i], (value >> i*8) & 0xff);
			}
		}
		
		static u_int64_t decodeNumeric(RangeDecoder& rangeDecoder, Order0Model& byteCountModel, vector<Order0Model>& numericModels){
			u_int8_t byteCount = rangeDecoder.nextByte(byteCountModel);
			//cout << (int)byteCount << endl;
			u_int64_t value = 0;
			for(int i=0; i<byteCount; i++){
				u_int8_t byteValue = rangeDecoder.nextByte(numericModels[i]);
				//cout << "Utils: " << (byteValue << i*8) << endl;
				value |= (byteValue << i*8);
			}
			
			return value;
		}
		
		static u_int64_t decodeFixedNumeric(RangeDecoder& rangeDecoder, vector<Order0Model>& numericModels, int byteCount){
			
			u_int64_t value = 0;
			for(int i=0; i<byteCount; i++){
				
				u_int8_t byteValue = rangeDecoder.nextByte(numericModels[i]);
				value |= (byteValue << i*8);
			}
			
			return value;
		}
		
		/*
		static stringstream _convert;
		
		static string numberToString(u_int64_t number){
			_convert << number;
			string result(_convert.str());
			//cout << result << endl;
			return result;
		}*/
		
		/*
		static int encodeDeltaValue(u_int64_t value, u_int64_t prevValue){
			int deltaType;

			int valueByteCount = CompressionUtils::getByteCount(value);
			
			//#ifdef PRINT_DEBUG_ENCODER
			//	cout << "\t\t\tPrev value: " << prevValue << endl;
			//	cout << "\t\t\tField value: " << value << "    Byte: " << valueByteCount << endl;
			//#endif
			
			u_int64_t deltaValue = value - prevValue;
			int deltaByteCount = CompressionUtils::getByteCount(deltaValue);
			
			//#ifdef PRINT_DEBUG_ENCODER
			//	cout << "\t\t\tDelta value: " << deltaValue << "    Byte: " << deltaByteCount << endl;
			//#endif
			if(deltaValue > 0 && deltaByteCount <= valueByteCount){
				//_rangeEncoder.encode(_typeModel[_misIndex], FIELD_DELTA);
				value = deltaValue;
				//valueByteCount = deltaByteCount;
				return 1;
			}
			else{
			
				
				u_int64_t deltaValue2 = prevValue - value;
				int deltaByteCount2 = CompressionUtils::getByteCount(deltaValue2);
			
				if(deltaValue2 > 0 && deltaByteCount2 <= valueByteCount){
					//_rangeEncoder.encode(_typeModel[_misIndex], FIELD_DELTA_2);
					value = deltaValue2;
					//valueByteCount = deltaByteCount2;
					return 2;
				}
				else{
					//_rangeEncoder.encode(_typeModel[_misIndex], FIELD_NUMERIC);
					return 0;
				}
			}
			
			
		}*/
	
};

	

	

#endif /* _COMPRESSIONUTILS_HPP_ */

