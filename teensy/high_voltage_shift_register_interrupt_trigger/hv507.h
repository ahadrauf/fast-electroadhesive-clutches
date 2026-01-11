// Starter code taken from https://github.com/Ai92/HV507/tree/main and https://github.com/DostoJoe/Microchip-HV507-Arduino-Library/tree/main
#ifndef HV507_H_
#define HV507_H_

#include <stdint.h>
#include "Arduino.h"

class HV507
{
public:
  void init(int8_t _clockPin, int8_t latchPin, int8_t dataPin, int8_t polarityPin, int8_t blankingPin, int8_t directionPin);
  void loadBitToRegister(uint8_t dataBit);  // 1 = on, 0 = off
	void loadByteToRegister(uint8_t data);
	void loadArrayToRegister(uint16_t data1, uint16_t data2, uint16_t data3, uint16_t data4);
	void loadArrayDirectlyToOutput(uint16_t data1, uint16_t data2, uint16_t data3, uint16_t data4);
	void loadRegisterToOutput();
	void loadInverseRegisterToOutput();
	void allOutputsLow();
	void allOutputsHigh();
	void invertOutputs();

private:
  int8_t _clockPin;
  int8_t _latchPin;
  int8_t _dataPin;
  int8_t _polarityPin;
  int8_t _blankingPin;
  int8_t _directionPin;
  bool _invertedOutput;
  bool _isInitialized;
};

#endif