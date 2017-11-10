/*
 * solver_bndy_isFocusBndy.cpp
 *
 *  Created on: Feb 6, 2014
 *      Author: chuan
 */

#include "../misc/misc_interpl.h"
#include "solver_fastSOR.h"

bool CDelphiFastSOR::isFocusBndy(delphi_real *** phimap)
{
   char toblbl[20],label[10],title[60],botlbl[16];
   delphi_real fScale1,goff,gmid;
   SGrid<delphi_real> fgBoxCenter1,fgXYZ1,fgXYZ2,fgTemp;
   delphi_integer iGrid1,isgrid;
   int ix,iy,iz,iout;

#ifndef MCCE
   ifstream ifPhiiFileStream;
   ifPhiiFileStream.open(strPhiiFile.c_str());
   if (!ifPhiiFileStream.is_open())
   {
      throw CUnknownPhiiFile(strPhiiFile);
      return false;
   }

   cout << "\n";
   cout << "focussing boundary condition\n";
   cout << "read from file\n";
   cout << strPhiiFile << endl << endl;

   ifPhiiFileStream.read(toblbl,21);
   ifPhiiFileStream.read(label ,10);
   ifPhiiFileStream.read(title ,54);

   for (iz = 0; iz < iGrid; iz++)
   {
      for (iy = 0; iy < iGrid; iy++)
      {
         for (ix = 0; ix < iGrid; ix++)
         {
            ifPhiiFileStream.read(reinterpret_cast<char*>(&phimap[iz][iy][ix]),sizeof(delphi_real));
         }
      }
   }

   ifPhiiFileStream.read(botlbl,16);
   ifPhiiFileStream.read(reinterpret_cast<char*>(&fScale1),sizeof(fScale1));
   ifPhiiFileStream.read(reinterpret_cast<char*>(&fgBoxCenter1),sizeof(fgBoxCenter1));
   ifPhiiFileStream.read(reinterpret_cast<char*>(&iGrid1),sizeof(iGrid1));

   ifPhiiFileStream.close();
#else
   prgfPhiMap   = pmcce->phimap;
   fScale1      = pmcce->scale1;
   fgBoxCenter1 = pmcce->oldmid1;
   iGrid1       = pmcce->igrid1;
#endif

   if (fZero > abs(fScale1-fScale))
   {
      cout << "scales are the same.\n";
      cout << "therefore assuming this to be a continuence\n";
   }
   else
   {
      cout << "\n";
      cout << " focussing potential map:\n";
      cout << title << endl;
      cout << "original scale (grids/A)      : " << fScale1 << endl;
      cout << "object centre at (A) : " << fgBoxCenter1.nX << fgBoxCenter1.nY << fgBoxCenter1.nZ << endl;

      //----- check to see that new grid lies within old one that is going to provide bc's
      iout = 0; goff = (iGrid1+1.0)/2.0;

      for (iz = 0; iz < iGrid; iz += iGrid-1)
      {
         for (iy = 0; iy < iGrid; iy += iGrid-1)
         {
            for (ix = 0; ix < iGrid; ix += iGrid-1)
            {
               fgXYZ1.nX = (delphi_real)(ix+1); fgXYZ1.nY = (delphi_real)(iy+1); fgXYZ1.nZ = (delphi_real)(iz+1);

               //for each new grid corner, calculate old grid coords
               fgXYZ2 = (fgXYZ1-goff)/fScale+fgBoxCenter;
               fgTemp = (fgXYZ2-fgBoxCenter1)*fScale1+goff;

               if (optORLE<delphi_real>(fgTemp,1.0) || optORGE<delphi_real>(fgTemp,(delphi_real)iGrid1)) iout = 1;
            }
         }
      }

      if (0 != iout) throw COutsideFocus(fScale1,fgBoxCenter1,fScale,fgBoxCenter);

      /*
       * for each boundary point
       *    convert to delphi_real coordinates
       *    convert to old grid coordinates
       *    interpolate potential
       * note that can use same potential array for boundaries since old potentials at boundary are not used for new ones
       */

      //----- save new grid size, and set temporarily to 65
      isgrid = iGrid; iGrid  = iGrid1; gmid = (isgrid+1.0)/2.0;
      cout << "pulling boundary values out of old potential map...\n";

      for (iz = 1; iz < isgrid-1; iz++)
      {
         for (iy = 1; iy < isgrid-1; iy++)
         {
            for (ix = 0; ix < isgrid; ix += isgrid-1)
            {
               fgXYZ1.nX = (delphi_real)(ix+1); fgXYZ1.nY = (delphi_real)(iy+1); fgXYZ1.nZ = (delphi_real)(iz+1);
               fgXYZ2 = (fgXYZ1-gmid)/fScale+fgBoxCenter;
               fgTemp = (fgXYZ2-fgBoxCenter1)*fScale1+goff;

               //for each new grid side, calculate old grid coords find potential
               phimap[iz][iy][ix] = interpl(iGrid,phimap,fgTemp);
            }
         }
      }

      for (iz = 1; iz < isgrid-1; iz++)
      {
         for (iy = 0; iy < isgrid; iy += isgrid-1)
         {
            for (ix = 1; ix < isgrid-1; ix++)
            {
               fgXYZ1.nX = (delphi_real)(ix+1); fgXYZ1.nY = (delphi_real)(iy+1); fgXYZ1.nZ = (delphi_real)(iz+1);
               fgXYZ2 = (fgXYZ1-gmid)/fScale+fgBoxCenter;
               fgTemp = (fgXYZ2-fgBoxCenter1)*fScale1+goff;

               //for each new grid side, calculate old grid coords find potential
               phimap[iz][iy][ix] = interpl(iGrid,phimap,fgTemp);
            }
         }
      }

      for (iz = 1; iz < isgrid; iz += isgrid-1)
      {
         for (iy = 1; iy < isgrid-1; iy++)
         {
            for (ix = 1; ix < isgrid-1; ix++)
            {
               fgXYZ1.nX = (delphi_real)(ix+1); fgXYZ1.nY = (delphi_real)(iy+1); fgXYZ1.nZ = (delphi_real)(iz+1);
               fgXYZ2 = (fgXYZ1-gmid)/fScale+fgBoxCenter;
               fgTemp = (fgXYZ2-fgBoxCenter1)*fScale1+goff;

               //for each new grid side, calculate old grid coords find potential
               phimap[iz][iy][ix] = interpl(iGrid,phimap,fgTemp);
            }
         }
      }

      iGrid = isgrid;
   }

   return true;
}
