#include<iostream>
#include "./map_zdcPix.C"


void test_map(){
  int Zdc_n=56;
  int Ch_max=24;
  int side=0;
  for(int ich=0;ich<Ch_max;ich++){
    
    double xcord, ycord, r, phi;
    set_xy_ZdcPix_had(ich, xcord, ycord);
    if(side==1) xcord=-xcord;
    set_rphi_ZdcPix(xcord, ycord, r, phi);
    cout<<"Ch ID = "<<ich<<",\t xcord ="<<xcord<<",\t ycord = "<<ycord<<",\t r = "<<r<<",\t phi = "<<phi<<endl;
  }
}


