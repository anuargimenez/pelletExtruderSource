   // Thermal conductivity kappa [W/(m.K)] function of T, for POLYSTYRENE (PS), 
   // Working range: 297.95K - 507.95K
    
    int i_index;
    scalar kappa_T_PS;
    scalar dT=5;
    scalar T0=297.95;
    scalar Temp_PS;

   // check that T is within a physical range
    if (T>750 || T<273.15)
      {
          FatalErrorIn
          (
              " Attempt to use polymerK,"
              " out of temperature range "
          )   << " T_kappa="<<T<<"; "
              << abort(FatalError);
      }
   

   // Thermal conductivity kappa [W/(m.K)] , PS :
   // data table extracted from J.Kattinger et al. (2022), Stuttgart.

    scalar kappa_PS[43]=
       {
	0.161937, 0.164796601602618, 0.167656203205237, 0.170515804807855,
	0.173375406410473, 0.176235008013092, 0.17909460961571, 0.181954211218329,
	0.185206122236239, 0.188496076422548, 0.191786030608857, 0.195075984795166,
	0.198365938981475, 0.201655893167784, 0.225431646436691, 0.251811996342789,
	0.249443695660769, 0.242542440347929, 0.241154123293213, 0.240812765955295,
	0.243137030020155, 0.245884457621725, 0.249732012339436, 0.253604260086491,
	0.254904692540174, 0.256205124993857, 0.25750555744754, 0.258805989901224,
	0.260106422354907, 0.26140685480859, 0.262707287262273, 0.264007719715957,
	0.263844631863019, 0.263444482970408, 0.263044334077796, 0.262644185185185,
	0.262244036292574, 0.261843887399963, 0.261443738507352, 0.26104358961474,
	0.260643440722129, 0.260243291829518, 0.259843142936907
       };
        
    
   // linear interpolation to calculate kappa(T) W/(m.K)
   
   // To avoid segmentation errors:     
      i_index = int(floor(fabs((T-T0)/dT)));
      
           if ( i_index > 42)
     {
     	i_index = 42;
     }
      
      Temp_PS = T0+i_index*dT;
      kappa_T_PS = kappa_PS[i_index]
                    + (kappa_PS[i_index+1]-kappa_PS[i_index])*(T-Temp_PS)/dT;

   //to avoid temperature limit errors
      if ( T > T0 + (dT*42))
    {
        return kappa_PS[42];
    }
    else if ( T < T0)
    {
        return kappa_PS[0];
    }
    else
    {
        return kappa_T_PS;
    }

   // end of kappa version implemented for PS




