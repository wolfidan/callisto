//  Easy_Diseqc V1.2, Monstein ETH Zurich, 20.06.2018.
// After an idea of University of Glasgow, UK. Private communication Graham W.
// Original code: https://github.com/acrerd/Arduino-Diseqc-solar-tracker
// DiSEqC-documentation: https://www.eutelsat.com/files/PDF/DiSEqC-documentation.zip

#include <Time.h>
#include <math.h>
#include  <util/parity.h>
#include <Wire.h>
#include <SoftwareSerial.h>

char Version[]="Easy_diseq.ino, 2019-10-11/cm"; // Version

float MaxRange = 69.0; // default maximum angular deflection depends on rotor type

int datapin1 = 8; // 22kHz Signal Pin
int datapin2 = 9; // 22kHz Signal Pin

float angle1 = 1;  // Azimuth/Hour angle as floating point number 0.0 ....+/-69
float angle2 = 1;  // Elevation/Declination as floating point number 0.0 ....+/-69
char tempbuf[80];  // keeps the command temporary until CRLF
String buffer;

void setup() {
 Serial.begin(9600);
 pinMode(datapin1,OUTPUT); // tone output
 pinMode(datapin2,OUTPUT); // tone output
}

void write01(){                      // write a '0' bit toneburst
 for (int i=1; i<=22; i++){         // 1 ms of 22 kHz (22 cycles)
  digitalWrite(datapin1,HIGH);
  delayMicroseconds(19);
  digitalWrite(datapin1,LOW);
  delayMicroseconds(19);
 }
 delayMicroseconds(500);             // 0.5 ms of silence
}

void write11(){                      // write a '1' bit toneburst
 for (int i=1; i<=11; i++){         // 0.5 ms of 22 kHz (11 cycles)
  digitalWrite(datapin1,HIGH);
  delayMicroseconds(20);
  digitalWrite(datapin1,LOW);
  delayMicroseconds(20);
 }
 delayMicroseconds(1000);            // 1 ms of silence
}

void write02(){                      // write a '0' bit toneburst
 for (int i=1; i<=22; i++){         // 1 ms of 22 kHz (22 cycles)
  digitalWrite(datapin2,HIGH);
  delayMicroseconds(19);
  digitalWrite(datapin2,LOW);
  delayMicroseconds(19);
 }
 delayMicroseconds(500);             // 0.5 ms of silence
}

void write12(){                      // write a '1' bit toneburst
 for (int i=1; i<=11; i++){         // 0.5 ms of 22 kHz (11 cycles)
  digitalWrite(datapin2,HIGH);
  delayMicroseconds(20);
  digitalWrite(datapin2,LOW);
  delayMicroseconds(20);
 }
 delayMicroseconds(1000);            // 1 ms of silence
}

void write_parity1(byte x){
 if (parity_even_bit(x)) write01(); else write11();
}

// write out a byte (as a toneburst)
// high bit first (ie as if reading from the left)
void write_byte1(byte x){
  for (int j=7; j>=0; j--){
  if (x & (1<<j)) write11(); else write01();
   }
}

void write_parity2(byte x){
 if (parity_even_bit(x)) write02(); else write12();
}

// write out a byte (as a toneburst)
// high bit first (ie as if reading from the left)
void write_byte2(byte x){
  for (int j=7; j>=0; j--){
  if (x & (1<<j)) write12(); else write02();
   }
}

// write out a byte with parity attached (as a toneburst)
void write_byte_with_parity1(byte x){
  write_byte1(x);
  write_parity1(x);
}

void write_byte_with_parity2(byte x){
  write_byte2(x);
  write_parity2(x);
}

void goto_angle1(float a){
  float fa16;
  byte n1,n2,n3,n4,d1,d2;
  int a16;
    
  if (a<0) { n1=0xE0;} else {n1=0xD0;}
  a16 =  (int) (16.0*abs(a)+0.5); 
  n2 = (a16 & 0xF00)>>8;
  d2 = a16 & 0xFF;
  d1 = n1 | n2;
  // send the command to the positioner
  noInterrupts();
  write_byte_with_parity1(0xE0);
  write_byte_with_parity1(0x31);
  write_byte_with_parity1(0x6E);
  write_byte_with_parity1(d1);
  write_byte_with_parity1(d2);
  interrupts();
}

void goto_angle2(float a){
  float fa16;
  byte n1,n2,n3,n4,d1,d2;
  int a16;
  
  if (a<0) { n1=0xE0;} else {n1=0xD0;}
  a16 =  (int) (16.0*abs(a)+0.5); 
  n2 = (a16 & 0xF00)>>8;
  d2 = a16 & 0xFF;
  d1 = n1 | n2;
  // send the command to the positioner
  noInterrupts();
  write_byte_with_parity2(0xE0);
  write_byte_with_parity2(0x31);
  write_byte_with_parity2(0x6E);
  write_byte_with_parity2(d1);
  write_byte_with_parity2(d2);
  interrupts();
}


void loop() {
        while (Serial.available() > 0)
        {
          int tmp;
          char st[20];
          char rx = Serial.read();  // read a single charecter
          buffer += rx;  //add character to the string buffer
          //Serial.print(rx);
         
          if ((rx == '\n') || (rx == '\r'))
          {
            buffer.toCharArray(tempbuf, 40);
            if (buffer.startsWith("max"))
            {
                sscanf(tempbuf,"max%s",&st); // extract attenuation as floating point number
                MaxRange = strtod(st,NULL);
                //Serial.println(MaxRange);
            } else
            if (buffer.startsWith("azi"))
            {
                sscanf(tempbuf,"azi%s",&st); // extract attenuation as floating point number
                angle1 = strtod(st,NULL);
                if (angle1 <-MaxRange) {angle1=-MaxRange;}
                if (angle1 > MaxRange) {angle1= MaxRange;}
                //Serial.println(angle1);
                goto_angle1(angle1);
            } else
            if (buffer.startsWith("ele"))
            {
                sscanf(tempbuf,"ele%s",&st); // extract attenuation as floating point number
                angle2 = strtod(st,NULL);
                if (angle2 <-MaxRange) {angle2=-MaxRange;}
                if (angle2 > MaxRange) {angle2= MaxRange;}
                //Serial.println(angle2);
                goto_angle2(angle2);
            } else      
            if (buffer.startsWith("-v")) // Check version of Arduino software
            {
               Serial.println(Version);
            } else
            if (buffer.startsWith("-h"))
            {
                Serial.println("Controller for two DISEqC satellite rotors.");
                Serial.println(Version);
                Serial.println("Command for max angle: maxDDD.d (floating point +/-°)");
                Serial.println("Command for elevation: eleDDD.d (floating point +/-°)");
                Serial.println("Command for azimuth  : aziDDD.d (floating point +/-°)");
                Serial.println("Extrem angles will be limited to MaxRange, stored in this program.");
                Serial.print  ("MaxRange=");
                Serial.println(MaxRange);
                Serial.println("-v sends software version");
                Serial.println("-h sends this help text");
            } 
            buffer = "";  //erase buffer for next command
          }
          
        }


}
