#include<iostream>

void set_grid_ZdcPix_had(const int chId, int& nX, int& nY) {

   nX = chId%4; nY = 5 - chId/4;
}

void set_xy_ZdcPix_had(const int chId, double& xcord, double& ycord) {

   int nx = chId%4;
   int ny = chId/4;
   
   if (ny == 0) {
      ycord = 4.0; 
      if      (nx == 0) xcord = -3.0;
      else if (nx == 1) xcord = -1.0;
      else if (nx == 2) xcord =  1.0;
      else if (nx == 3) xcord =  3.0;
   }

   else if (ny == 1) {
      ycord = 6.5/3.0;
      if      (nx == 0) xcord = -8.5/3.0;
      else if (nx == 1) xcord = -2.5/3.0;
      else if (nx == 2) xcord =  3.5/3.0;
      else if (nx == 3) xcord =  9.5/3.0;
   }

   else if (ny == 2) {
      ycord = 2.5/3.0;
      if      (nx == 0) xcord = -9.5/3.0;
      else if (nx == 1) xcord = -3.5/3.0;
      else if (nx == 2) xcord =  2.5/3.0;
      else if (nx == 3) xcord =  8.5/3.0;
   }

   else if (ny == 3) {
      ycord = -2.5/3.0;
      if      (nx == 0) xcord = -8.5/3.0;
      else if (nx == 1) xcord = -2.5/3.0;
      else if (nx == 2) xcord =  3.5/3.0;
      else if (nx == 3) xcord =  9.5/3.0;
   }
   
   else if (ny == 4) {
      ycord = -6.5/3.0;
      if      (nx == 0) xcord = -9.5/3.0;
      else if (nx == 1) xcord = -3.5/3.0;
      else if (nx == 2) xcord =  2.5/3.0;
      else if (nx == 3) xcord =  8.5/3.0;
   }
   
   else if (ny == 5) {
      ycord = -4.0;
      if      (nx == 0) xcord = -3.0;
      else if (nx == 1) xcord = -1.0;
      else if (nx == 2) xcord =  1.0;
      else if (nx == 3) xcord =  3.0;
   }

   else {
      ycord = 1000.0;
      xcord = 1000.0;
      std::cout << "Unknown Channel in set_xy_ZdcPix_had!! chId = " << chId << std::endl;
   }
}

void get_gridFull_ZdcPix_had(const int chId, std::vector<int>& gridX, std::vector<int>& gridY) {

   gridX.assign(4,-1);
   gridY.assign(4,-1);
   static int gridMap[24][4] = {{0 ,1 ,8 ,9 }, {2 ,3 ,10,11}, {4 ,5 ,12,13}, {6 ,7 ,14,15},
                                {16,17,25,-1}, {18,19,27,-1}, {20,21,29,-1}, {22,23,31,-1},
                                {24,32,33,-1}, {26,34,35,-1}, {28,36,37,-1}, {30,38,39,-1},
                                {40,41,49,-1}, {42,43,51,-1}, {44,45,53,-1}, {46,47,55,-1},
                                {48,56,57,-1}, {50,58,59,-1}, {52,60,61,-1}, {54,62,63,-1},
                                {64,65,72,73}, {66,67,74,75}, {68,69,76,77}, {70,71,78,79}
                               }; 

   for (int ig=0; ig<4; ig++) {
      if (gridMap[chId][ig]<0) continue;
      int x = (gridMap[chId][ig])%8;
      int y = (gridMap[chId][ig])/8;
      gridX.at(ig) = x;
      gridY.at(ig) = 9-y;
   }
}

void set_rphi_ZdcPix_had (const int chId, double& r, double& phi) {

   double xcord, ycord;
   set_xy_ZdcPix_had(chId, xcord, ycord);

   if (fabs(xcord-1000.0) < 1e-6 && fabs(ycord-1000.0) < 1e-6) {
      r   = 1000.0;
      phi = 1000.0;
      //std::cout << "Unknown Channel in set_rphi_ZdcPix_had!! chId = " << chId << std::endl;
      return;
   }

   r = sqrt( xcord*xcord + ycord*ycord );
   phi = atan2(ycord, xcord);
}

void set_rphi_ZdcPix (const double xcord, const double ycord, double& r, double& phi) {

   r = sqrt( xcord*xcord + ycord*ycord );
   phi = atan2(ycord, xcord);
}


void set_grid_ZdcPix_em(const int chId, int& nX, int& nY) {

      nX = chId%8; nY = 7 - chId/8;
}

void set_xy_ZdcPix_em(const int chId, double& xcord, double& ycord) {

   int nx = chId%8;
   int ny = chId/8;

   if (nx == 0)      xcord = -3.5;
   else if (nx == 1) xcord = -2.5;
   else if (nx == 2) xcord = -1.5;
   else if (nx == 3) xcord = -0.5;
   else if (nx == 4) xcord =  0.5;
   else if (nx == 5) xcord =  1.5;
   else if (nx == 6) xcord =  2.5;
   else if (nx == 7) xcord =  3.5;
  
   if (ny == 0)      ycord =  3.5;
   else if (ny == 1) ycord =  2.5;
   else if (ny == 2) ycord =  1.5;
   else if (ny == 3) ycord =  0.5;
   else if (ny == 4) ycord = -0.5;
   else if (ny == 5) ycord = -1.5;
   else if (ny == 6) ycord = -2.5;
   else if (ny == 7) ycord = -3.5; 
}
