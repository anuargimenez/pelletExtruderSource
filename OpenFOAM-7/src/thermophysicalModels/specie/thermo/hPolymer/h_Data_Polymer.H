    // enthalpy [J/kg] function of T, for PS,
    // tabulated for T from T0=293.362K to 513.362K 
    // with tabulation interval of dT=5K

    int i_index;
    scalar dT=5;
    scalar T0=293.362;
    scalar Temp_PS;
    scalar h_T_PS;

     // check that T is within a physical range
     if ( T>750 || T<273.15)
        {
            FatalErrorIn
            (
              " Attempt to use hPolymerThermo,"
              " out of temperature range "
            )   << " T_h="<<T<<"; "
                << abort(FatalError);
       }
     
     // enthalpy [J/kg] , PS : 
     // data table from J.Kattinger et al. (2022), Stuttgart.
     // so that dh/dT=Cp as needed for calculating T
     
     scalar h_PS[45] = 
     {
    0, 6116.50021374658, 12474.6008549863, 19074.2689041986, 25909.4361749525, 32925.0787670808, 40072.2079473265, 47352.3522493183, 54769.8894816598, 62341.76243174, 70082.0646119726, 77990.7960223575, 86065.3165970228, 94284.9685435426, 102631.734135363, 111221.12235733, 120371.482893105, 130000.568886676, 139725.412043795, 149508.737654016, 159408.980636124, 169401.299912449, 179463.025515511, 189595.487409517, 199802.865243136, 210088.008700832, 220450.917782603, 230890.472667892, 241405.476455715, 251997.141890989, 262671.736346596, 273434.23737008, 284284.237555923, 295201.959937263, 306168.034952758, 317182.462602408, 328245.242886211, 339356.375804169, 350515.861356281, 361728.125271616, 372998.142515405, 384326.462323807, 395713.084696824, 407158.009634455, 418661.2371367
      };
         
     // linear interpolation to calculate h(T) in J/kg
     
     // to avoid segmentation errors
     i_index = int(floor(fabs((T-T0)/dT)));
     
     if ( i_index > 44)
     {
     	i_index = 44;
     }
     
     Temp_PS = T0+i_index*dT;
     h_T_PS = h_PS[i_index]
               + (h_PS[i_index+1]-h_PS[i_index])*(T-Temp_PS)/dT;

    // to avoid limit errors
    if ( T > T0 + (dT*44))
    {
        return h_PS[44];
    }
    else if ( T < T0)
    {
        return h_PS[0];
    }
    else
    {
        return h_T_PS;
    }

   // end of h version implemented for PS 
              
