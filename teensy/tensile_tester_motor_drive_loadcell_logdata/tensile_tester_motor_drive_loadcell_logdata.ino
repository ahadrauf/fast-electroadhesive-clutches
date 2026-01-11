/***
Teensy 4.1 code for the tensile testing setup

The Teensy microcontroller pre-programs the desired max motor acceleration to a Tic 36v4 H-Bridge motor 
    driver (Pololu).
When it receives a Serial trigger over USB (two numbers: the distance to move, and then the speed), it starts
    driving the stepper motor at a constant speed forwards. It also triggers the load cell ADC and saves
    the load cell reading and current high voltage reading into an array. The array is saved onto additional
    PSRAM soldered onto the Teensy 4.1 (https://www.pjrc.com/store/psram.html). After the tensile testing run
    is finished, the Teensy uploads all the data in bulk back to the host PC.
***/

#include <SPI.h>
#include <Tic.h>
#include "ADS1263.h"

// Motor driver I2C pins
const int DIR_PIN = 37;
const int STEP_PIN = 36;
const int CS_PIN = 10;
const int INTERRUPT_PIN = 22;
const int VOLTAGE_IN_PIN = A2;
const int HV_ENABLED_INTERRUPT_PIN = 32;
TicI2C tic;

// Motor driver settings
const uint32_t NEMA14_STEPS_PER_REV = 200;
const uint32_t NEMA14_LEAD_SCREW_PITCH = 2;  // mm
const uint32_t NEMA14_MICROSTEPS = 256;
const uint32_t NEMA14_STEPS_PER_MM = NEMA14_STEPS_PER_REV * NEMA14_MICROSTEPS / NEMA14_LEAD_SCREW_PITCH; // 25600
const uint32_t MMS_TO_TIC_SPEED = NEMA14_STEPS_PER_MM * 10000;  // pulses/10000 s
const uint32_t MMS2_TO_TIC_ACCEL = NEMA14_STEPS_PER_MM * 100;  // pulses/s / 100s
const float maxSpeed = 0.5;  // mm/s
const float maxAccel = 10;  // mm/s^2
int32_t targetPosition = 0;
uint8_t LoadCellGain = 32;
const float LC2F_SLOPE = 64.69523324;
const float LC2F_BIAS = -0.112957877;

// Stepper motor linear tensile testing parameters
float distance_first_movement;
float time_first_movement = 2;  // sec
float distance_brake_engaged;
float speed_brake_engaged;
bool engage_brake_early;
unsigned int delay_before_engaging_brake = 1000;  // ms
unsigned int delay_while_engaging_brake = 1000;  // ms
unsigned int delay_after_releasing_brake = 2000;  // ms
float distance_last_movement;
float time_brake_released;
float time_end_forward_movement;
float time_end_backward_movement;
float distance_total;
int motion_state = 0;  // 0 = not in motion, 1 = moving forward, 2 = waiting, 3 = moving back
bool brake_engaged = false;
elapsedMillis time_engaging_brake;
elapsedMillis time_before_engaging_brake;
const float distance_multiplier = 0.1;
const float speed_multiplier = 0.05;
volatile bool hv_enabled = false;

// Data logging settings
const uint32_t SerialUSBSpeed = 250000;
const uint8_t csPin = 10;
const uint8_t dataReadyPin = 9;
volatile double adcReading = 0.0;
volatile bool newADCReading = false;
ADS1263 adc(csPin, dataReadyPin);
bool record_data = false;
elapsedMicros time_since_run_start;
EXTMEM unsigned long times[250000];
EXTMEM float loadcell_forces[250000];
EXTMEM float moving_forwards[250000];
EXTMEM float voltage_out[250000];
unsigned long currIdx = 0;

void setup() {
  Serial.begin(115200);
  Serial.setTimeout(1);

  pinMode(VOLTAGE_IN_PIN, INPUT);
  pinMode(INTERRUPT_PIN, OUTPUT);
  pinMode(HV_ENABLED_INTERRUPT_PIN, INPUT);
  digitalWrite(INTERRUPT_PIN, LOW);
  currIdx = 0;
  motion_state = 0;

  // Initializes the motor driver
  Wire.begin();
  delay(20);
  tic.haltAndSetPosition(0); delay(20);
  tic.setMaxAccel((uint32_t) (maxAccel * MMS2_TO_TIC_ACCEL)); delay(20);
  tic.setMaxDecel((uint32_t) (maxAccel * MMS2_TO_TIC_ACCEL)); delay(20);
  tic.enterSafeStart(); delay(20);

  // Initializes the ADS1263 device
  SPI.setSCK(13);
  SPI.begin();
  adc.BeginCustom(LoadCellGain);
  adc.StartADC1();
  attachInterrupt(digitalPinToInterrupt(dataReadyPin), ADC1Callback, RISING);
  attachInterrupt(digitalPinToInterrupt(HV_ENABLED_INTERRUPT_PIN), setHVEnabled, CHANGE);
}

void loop() {
  if (motion_state == 0 && Serial.available() >= 2) {
    distance_brake_engaged = distance_multiplier * Serial.read(); // distance to move while braked
    speed_brake_engaged = speed_multiplier * Serial.read(); // speed to move at while braked
    distance_last_movement = delay_after_releasing_brake * speed_brake_engaged / 1000.0;
    
    record_data = false;
    motion_state = 1;
    time_engaging_brake = 0;
    digitalWrite(INTERRUPT_PIN, LOW);
    resetArrays();
    tic.exitSafeStart(); delay(100);
    tic.setMaxSpeed((uint32_t) (speed_brake_engaged * MMS_TO_TIC_SPEED)); delay(20);
    tic.setStartingSpeed((uint32_t) (speed_brake_engaged * MMS_TO_TIC_SPEED)); delay(20);

    distance_first_movement = time_first_movement * speed_brake_engaged;  // distance_last_movement
    distance_total = distance_first_movement + distance_brake_engaged;
    time_brake_released = time_first_movement + distance_brake_engaged / speed_brake_engaged;
    time_end_forward_movement = time_brake_released + delay_after_releasing_brake / 1000;
    time_end_backward_movement = time_end_forward_movement + distance_total / speed_brake_engaged;
    record_data = true;
    time_since_run_start = 0;
    setTargetPosition(distance_total);
  }

  if (motion_state == 1 && (time_since_run_start >= time_first_movement * 1e6)) {
    motion_state = 4;
    time_engaging_brake = 0;
    brake_engaged = true; 
    digitalWrite(INTERRUPT_PIN, HIGH);
  } else if (motion_state == 4 && (time_since_run_start >= time_brake_released * 1e6)) {
    motion_state = 6;  // start moving forward on final movement
    brake_engaged = false;
    digitalWrite(INTERRUPT_PIN, LOW);
    digitalWrite(INTERRUPT_PIN, LOW);
    time_engaging_brake = 0;
  } else if (motion_state == 6 && (time_since_run_start >= time_end_forward_movement * 1e6)) {  // after moving while the brake is engaged
    motion_state = 7;
    record_data = false;
    setTargetPosition(0);
  } else if (motion_state == 7 && (time_since_run_start >= time_end_backward_movement * 1e6)) {  // after moving while the brake is engaged
    motion_state = 0;
    tic.enterSafeStart(); delay(20);
    cli();
    writeDataToSerial();
    sei();
  }
  
  if (record_data && newADCReading) { 
    times[currIdx] = time_since_run_start;
    loadcell_forces[currIdx] = LC2F_SLOPE*adcReading + LC2F_BIAS;
    newADCReading = false;
    moving_forwards[currIdx] = hv_enabled;
    voltage_out[currIdx] = 0.344*analogRead(VOLTAGE_IN_PIN);
    currIdx++;
  }
}

void setTargetPosition(float position_mm) {
  targetPosition = (int32_t) (position_mm * NEMA14_STEPS_PER_MM);
  tic.setTargetPosition(targetPosition);
}

void ADC1Callback() {
  adc.ReadADC1();
	adcReading = adc.GetADC1ValueAfterPGA();
  newADCReading = true;
}

void writeDataToSerial(void) {
  Serial.println(currIdx);
  for (unsigned long i = 0; i < currIdx; i++) { Serial.println(times[i]); delayMicroseconds(50); }
  Serial.println(currIdx);
  for (unsigned long i = 0; i < currIdx; i++) { Serial.println(loadcell_forces[i], 8); delayMicroseconds(50); }
  Serial.println(currIdx);
  for (unsigned long i = 0; i < currIdx; i++) { Serial.println(moving_forwards[i]); delayMicroseconds(50); }
  Serial.println(currIdx);
  for (unsigned long i = 0; i < currIdx; i++) { Serial.println(voltage_out[i], 6); delayMicroseconds(50); }
  Serial.println(currIdx);
}

void resetArrays() {
  currIdx = 0;
  for (int i = 0; i < 250000; i++) {
    times[i] = 0;
    loadcell_forces[i] = 0;
    moving_forwards[i] = 0;
    voltage_out[i] = 0;
  }
}

void setHVEnabled() {
  hv_enabled = digitalRead(HV_ENABLED_INTERRUPT_PIN);
}