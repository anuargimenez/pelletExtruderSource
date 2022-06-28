// Density [kg/m3] (specific volume m3/kg) function of T and P, for PS

// TAIT EQ. OF STATE for amorphous polymers:
    scalar spec_volS;
    scalar spec_volM;

    //melt
     spec_volM = (b1m_+b2m_*((-b6_*p)-b5_+T))*(1-0.0894*log(1+p*pow(b3m_,-1)*exp(b4m_*((-b6_*p)-b5_+T))));
     
    //solid
     spec_volS = (b1s_+b2s_*((-b6_*p)-b5_+T))*(1-0.0894*log(1+p*pow(b3s_,-1)*exp(b4s_*((-b6_*p)-b5_+T))));

  
// INTERPOLATED TRANSITION REGION
    scalar Td = 10.0; //[K] half the interpolation band
    scalar Tt = b5_ + b6_*p; //Tait solidification(transition) temperature
    scalar Tp = Tt + Td; //upper band limit
    scalar Tm = Tt - Td; //lower band limit

    if ( T > Tp ) // melt
    {
        return 1.0/spec_volM;
    }
    else if ( T > Tm) // interpolated transition region
    {
        return 1.0/((T-Tm)/(Tp-Tm) * spec_volM + (Tp-T)/(Tp-Tm) * spec_volS);
    }
    else // solid
    {
        return 1.0/spec_volS;
    }

