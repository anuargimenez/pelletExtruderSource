    // heat capacity at constant pressure [J/(kg.K)] function of T, for PS
    // tabulated for T from T0=293.362K to 513.362K 
    // with tabulation interval of dT=5K

    int i_index;
    scalar dT=5;
    scalar T0=293.362;
    scalar Temp_PS;
    scalar Cp_T_PS;

    // check that T is within a physical range
    if (T>750 || T<273.15)
        {
            FatalErrorIn
            (
                " Attempt to use hPolymerThermo,"
                " out of temperature range "
            )   << " T_Cp="<<T<<"; "
                << abort(FatalError);
        }
               

     // Specific heat at constant pressure (J/(kg K)), PS : 
     // data table extracted from J.Kattinger et al. (2022), Stuttgart.
     // so that dh/dT=Cp as needed for calculating T
     
      scalar Cp_PS[45]=
   {
   1199.14, 1247.46008549863, 1295.78017099726, 1344.08704868764, 1389.97985961394, 1416.2771772374, 1442.57449486087, 1469.48322593583, 1497.5316670008, 1531.21751303128, 1564.90335906175, 1598.58920509222, 1631.21902477389, 1656.64175383405, 1682.06448289422, 1753.69080589255, 1906.45340841731, 1945.18098901099, 1944.75627383677, 1968.57397025172, 1991.52322259136, 2005.40448793876, 2019.28575328615, 2033.69900431621, 2049.25212913144, 2064.80525394667, 2080.3583787619, 2095.46357535351, 2110.53793977581, 2126.12823433379, 2143.70954790908, 2161.29086148437, 2178.70921285272, 2188.37973968358, 2198.05026651444, 2207.72079334529, 2217.39132017615, 2227.06184700701, 2236.73237383787, 2248.17319229618, 2259.83370521904, 2271.4942181419, 2283.15473106476, 2294.81524398763, 2306.47575691049
     };
         
     // linear interpolation to calculate Cp(T) in J/(kg.K)
     
     // to avoid segmentation errors     
     i_index = int(floor(fabs((T-T0)/dT)));
     
          if ( i_index > 44)
     {
     	i_index = 44;
     }
     
     Temp_PS = T0+i_index*dT;
     Cp_T_PS = Cp_PS[i_index]
                + (Cp_PS[i_index+1]-Cp_PS[i_index])*(T-Temp_PS)/dT;

     
    // to avoid limit errors
    if ( T > T0 + (dT*44))
    {
        return Cp_PS[44];
    }
    else if ( T < T0)
    {
        return Cp_PS[0];
    }
    else
    {
        return Cp_T_PS;
    }


   // end of cp version implemented for PS 
