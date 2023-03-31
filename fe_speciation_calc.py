import pandas as pd
import numpy as np
import periodictable as pt





def normalize(df, comps):
    """
    Normalize oxides data (weight percentage or moles)

    Parameters
    -----------

    df : :class:`pandas.DataFrame`
        dataframe to be normalized

    comps : :class:`list`
        list of optional oxides to be normalized
    Returns
    -------
        dfmc : :class:`pandas.DataFrame`
        normalized dataframe
    """
    dfmc = df.copy(deep = True)
    if comps:
        cmpnts = [c for c in comps if c in dfmc.columns]
        dfmc.loc[:, cmpnts] = 100 * dfmc.loc[:, cmpnts].divide(
            dfmc.loc[:, cmpnts].sum(axis=1).replace(0, np.nan), axis=0
            )
    else:
        dfmc = 100 * dfmc.divide(dfmc.sum(axis=1).replace(0, 100.0), axis=0)
    return dfmc





def weight2mole(df, comps=None, norm_fractions = True):
    """
    Transform weight percentage to mole (fractions)

    Parameters
    -----------

        df : :class:`pandas.DataFrame`
            dataframe in weight percentage to convert to mole fractions

        comps : :class:`list`
            list of optional oxides to be converted

        norm_fractions : :class:`boolean`
            default is True
            if norm_fractions:
                sum and convert moles to 100
            else:
                moles equal division of oxides and total.
    Returns
    -------
        :class:`pandas.DataFrame`
        converted dataframe
    """
    if comps == None:
            comps = ['SiO2', 'TiO2', 'Al2O3', 'FeO', 'MnO', 'MgO', 'CaO', 'P2O5', 'Na2O', 'K2O']
    dfmc = df.copy()
    dfmc = dfmc[comps]
    owms = [pt.formula(c).mass for c in dfmc.columns]
    if norm_fractions:
        return normalize(dfmc.div(owms), comps)
    else:
        return dfmc.div(owms)





def weight2mole_sc(df, comps=None, norm_fractions = True):
    """
    Transform weight percentage to mole (fractions) calculated on the single cation basis

    Parameters
    -----------

        df : :class:`pandas.DataFrame`
            dataframe in weight percentage to convert to mole fractions

        comps : :class:`list`
            list of optional oxides to be converted

        norm_fractions : :class:`boolean`
            default is True
            if norm_fractions:
                sum and convert moles to 100
            else:
                moles equal division of oxides and total.
    Returns
    -------
        :class:`pandas.DataFrame`
        converted dataframe
    """
    if comps == None:
            comps = ['SiO2', 'TiO2', 'Al2O3', 'FeO', 'MnO', 'MgO', 'CaO', 'Na2O', 'K2O', 'P2O5', 'Cr2O3', 'NiO']
    dfmc = df.copy()
    dfmc = dfmc[comps]
    dfmc[['Al2O3', 'Na2O', 'K2O', 'P2O5', 'Cr2O3']] = 2*dfmc[['Al2O3', 'Na2O', 'K2O', 'P2O5', 'Cr2O3']]
    ox_mass = [pt.formula(c).mass for c in dfmc.columns]
    if norm_fractions:
        return normalize(dfmc.div(ox_mass), comps)
    else:
        return dfmc.div(ox_mass)
    




def fe_speciation(liq_comps, index, T_K, P_Pa, fO2, model='Kress_1991', df_normalize=True):
    """
    Calculate Fe3+Fetot, FeO and Fe2O3 from liquid comps and fO2

    Parameters
    ----------
        liq_comps: pandas.DataFrame
            Liquid compositions in wt% with headings "SiO2", "MgO" etc
        
        index: pandas.Series
            Sample index
        
        T_K: int, flt, pandas.Series
            Temperature in Kelvin

        P_Pa: int, flt, pandas.Series
            Pressure in pascals
        
        fO2: int, flt, pandas.Series
            fO2 in bars

        model: str
            "IrvineBarager_1971" - used for whole rock data
            "LeMaitre_1976"
            "Sack_1980"
            "Kress_1991" - eqn7
            "Kress_1991_lowp" - eqn6
            "Jayasuriya_2004"
            "Putirka_2016" - eqn 6b
            "ONeill_2018" - eqn 9a
            "Borisov_2018"

        normalize: bool
            Default = True
            if true:
                Normalizes compositions to 100
            else:
                Returns unormalized dataframe with new total
        
    Returns
    -------
    df_out: 
            Liquid compositions with new FeO and Fe2O3 values, XFe2O3/XFeO and Fe3+/Fetot
    """

    #Comps
    comps = ['SiO2', 'TiO2', 'Al2O3', 'FeO', 'MnO', 'MgO', 'CaO', 'Na2O', 'K2O', 'P2O5', 'Cr2O3', 'NiO']
    comps_noFe = ['SiO2', 'TiO2', 'Al2O3', 'MnO', 'MgO', 'CaO', 'Na2O', 'K2O', 'P2O5', 'Cr2O3', 'NiO']

    comps_borisov = ['SiO2', 'TiO2', 'Al2O3', 'FeO', 'MgO', 'CaO', 'Na2O', 'K2O', 'P2O5']

    comps_normalize = ['SiO2', 'TiO2', 'Al2O3', 'FeO', 'Fe2O3', 'MnO', 'MgO', 'CaO', 'Na2O', 'K2O', 'P2O5', 'Cr2O3', 'NiO']
    comps_final = ['SiO2', 'TiO2', 'Al2O3', 'FeO_calc', 'Fe2O3_calc', 'MnO', 'MgO', 'CaO', 'Na2O', 'K2O', 'P2O5', 'Cr2O3', 'NiO']

    #Calculate normalised wt% df, mole fractions and mole fractions on single cation basis
    df = liq_comps.copy()
    df_moles = weight2mole(df, comps=comps, norm_fractions = True)/100
    df_moles_borisov = weight2mole(df, comps=comps_borisov, norm_fractions = True)/100 
    df_moles_sc = weight2mole_sc(df, comps=comps, norm_fractions = True)/100
    df_norm = normalize(df, comps)

    #Create empty output dataframe
    df_out = pd.DataFrame()



    if model == "IrvineBarager_1971":
        
        df_norm.loc[df_norm['FeO']*1.11111 > df_norm['TiO2']+1.5, 'Fe2O3_calc'] = df_norm['TiO2'] + 1.5
        df_norm['FeO_calc'] = (df_norm['FeO']*1.11111 - df_norm['Fe2O3_calc'])/1.11111

        if df_normalize == False:
            df_out = df[comps_noFe].copy()
            df_out['Fe2O3_calc'] = df_norm['Fe2O3_calc'].copy()
            df_out['FeO_calc'] = df_norm['FeO_calc'].copy()
            df_out = df_out[comps_final]
            df_out['Old_Total'] = pd.DataFrame.sum(df[comps], axis=1)
            df_out['New_Total'] = pd.DataFrame.sum(df[comps_noFe], axis=1) + df_out['Fe2O3_calc'] + df_out['FeO_calc']
            df_out['Fe3_Fetot'] = 1/((df_out.FeO_calc/(df_out.Fe2O3_calc*0.8997))+1)

        if df_normalize == True:
            df_out = (normalize(df[comps_noFe].copy(), comps_noFe)).multiply(((100-df_norm['Fe2O3_calc']-df_norm['FeO_calc'])/100), axis=0)
            df_out['Fe2O3_calc'] = df_norm['Fe2O3_calc'].copy()
            df_out['FeO_calc'] = df_norm['FeO_calc'].copy()
            df_out = df_out[comps_final]
            df_out['total'] = df_out.sum(axis=1)
            df_out['Fe3Fetot'] = 1/((df_out.FeO_calc/(df_out.Fe2O3_calc*0.8997))+1)





    elif model == "LeMaitre_1976":
        FeO_FeOFe2O3 = (0.93 - 0.0042*df_norm.SiO2 
                        - 0.022*(df_norm.Na2O + df_norm.K2O)
                        )
        Fe3_Fe2 = ((1/FeO_FeOFe2O3)-1)*2*71.8444/159.6882
        Fe3_Fetot = 1/((1/Fe3_Fe2)+1)
        Fe2O3_FeO = Fe3_Fe2*(159.69/(2*71.844))

        df_out['Fe2O3_calc'] = (Fe2O3_FeO*df.FeO/(0.8998*Fe2O3_FeO + 1))
        df_out['FeO_calc'] = df.FeO - 0.8998*df_out.Fe2O3_calc

        if df_normalize == False:
            df_out = df[comps_noFe].copy()
            df_out['Fe2O3_calc'] = (Fe2O3_FeO*df.FeO/(0.8998*Fe2O3_FeO + 1))
            df_out['FeO_calc'] = df.FeO - 0.8998*df_out.Fe2O3_calc
            df_out = df_out[comps_final]
            df_out['Old_Total'] = pd.DataFrame.sum(df[comps], axis=1)
            df_out['New_Total'] = pd.DataFrame.sum(df[comps_noFe], axis=1) + df_out['Fe2O3_calc'] + df_out['FeO_calc']
            df_out['Fe3Fetot'] = Fe3_Fetot

        if df_normalize == True:
            df_out['Fe2O3_calc'] = (Fe2O3_FeO*df.FeO/(0.8998*Fe2O3_FeO + 1))
            df_out['FeO_calc'] = df.FeO - 0.8998*df_out.Fe2O3_calc
            df_out[comps_noFe] = (normalize(df[comps_noFe].copy(), comps_noFe)).multiply(((100-df_out['Fe2O3_calc']-df_out['FeO_calc'])/100), axis=0)
            df_out = df_out[comps_final]
            df_out['total'] = df_out.sum(axis=1)
            df_out['Fe3Fetot'] = Fe3_Fetot





    elif model == "Sack_1980":

        #define variables
        a = 0.218130
        b = 13184.7
        c = -4.49933
        d_sio2 = -2.15036
        d_al2o3 = -8.35163
        d_feo = -4.49508
        d_mgo = -5.43639
        d_cao = 0.073113
        d_na2o = 3.54148
        d_k2o = 4.18688

        #input into eqn
        XFe2O3_XFeO = np.exp(a*np.log(fO2)
            + np.divide(b,(T_K))
            + c 
            + d_sio2*df_moles.SiO2 
            + d_al2o3*df_moles.Al2O3 
            + d_feo*df_moles.FeO
            + d_mgo*df_moles.MgO
            + d_cao*df_moles.CaO
            + d_na2o*df_moles.Na2O
            + d_k2o*df_moles.K2O
            )
    
        XFe2O3 = XFe2O3_XFeO*df_moles.FeO/(2*XFe2O3_XFeO + 1)
        XFeO = df_moles['FeO'] - 2*XFe2O3

        sum_moles = (df_moles['SiO2']*60.084 + df_moles['TiO2']*79.866 + df_moles['Al2O3']*101.96 + df_moles['MnO']*70.9374 
                + df_moles['MgO']*40.3044 + df_moles['CaO']*56.0774 + df_moles['P2O5']*141.943 + df_moles['Na2O']*61.9789 
                + df_moles['K2O']*94.2 + df_moles['Cr2O3']*151.9904 + df_moles['NiO']*74.6928 + XFeO*71.844 + XFe2O3*159.69
                )

        #New Wt% values
        Fe2O3_wt = 100*(XFe2O3*159.69)/sum_moles
        FeO_wt = 100*(XFeO*71.844)/sum_moles



        if df_normalize == False:
            df_out = df[comps_noFe].copy()
            df_out['Fe2O3_calc'] = Fe2O3_wt
            df_out['FeO_calc'] = FeO_wt
            df_out = df_out[comps_final]
            df_out['Old_Total'] = pd.DataFrame.sum(df[comps], axis=1)
            df_out['New_Total'] = pd.DataFrame.sum(df[comps_noFe], axis=1) + Fe2O3_wt + FeO_wt
            df_out['XFe3XFe2'] = XFe2O3_XFeO
            df_out['Fe3Fetot'] = 2*XFe2O3/(2*XFe2O3+XFeO)

        if df_normalize == True:
            df_out = pd.DataFrame(data={'SiO2': 100*df_moles['SiO2']*60.084/sum_moles,
                                    'TiO2': 100*df_moles['TiO2']*79.866/sum_moles,
                                    'Al2O3':100*df_moles['Al2O3']*101.96/sum_moles,
                                    'Fe2O3_calc': Fe2O3_wt,
                                    'FeO_calc': FeO_wt,
                                    'MnO': 100*df_moles['MnO']*70.9374/sum_moles,
                                    'MgO': 100*df_moles['MgO']*40.3044/sum_moles,
                                    'CaO': 100*df_moles['CaO']*56.0774/sum_moles,
                                    'Na2O': 100*df_moles['Na2O']*61.9789/sum_moles,
                                    'K2O': 100*df_moles['K2O']*94.2/sum_moles,
                                    'P2O5':  100*df_moles['P2O5']*141.943/sum_moles,
                                    'Cr2O3': 100*df_moles['Cr2O3']*151.9904/sum_moles,
                                    'NiO': 100*df_moles['NiO']*74.6928/sum_moles
                                    })

            df_out = df_out[comps_final]
            df_out['total'] = df_out.sum(axis=1)
            df_out['Fe3Fetot'] = 2*XFe2O3/(2*XFe2O3+XFeO)
            df_out['XFe2O3_XFeO'] = XFe2O3_XFeO





    elif model == "Kress_1991":
        T_o = 1673

        XFe2O3_XFeO = np.exp(0.196*np.log(fO2)
                + 11492/T_K
                - 6.675
                - 2.243*df_moles['Al2O3']
                - 1.828*df_moles['FeO']
                + 3.201*df_moles['CaO']
                + 5.854*df_moles['Na2O']
                + 6.215*df_moles['K2O']
                - 3.36*(1-T_o/T_K-np.log(T_K/T_o))
                - 0.000000701*P_Pa/T_K
                - 0.000000000154*(T_K-T_o)*P_Pa/T_K
                + 0.0000000000000000385*P_Pa**2/T_K
                )
        
        #Calculate new mole fractions
        XFe2O3 = XFe2O3_XFeO*df_moles.FeO/(2*XFe2O3_XFeO + 1)
        XFeO = df_moles['FeO'] - 2*XFe2O3


        sum_moles = (df_moles['SiO2']*60.084 + df_moles['TiO2']*79.866 + df_moles['Al2O3']*101.96 + df_moles['MnO']*70.9374 
                + df_moles['MgO']*40.3044 + df_moles['CaO']*56.0774 + df_moles['P2O5']*141.943 + df_moles['Na2O']*61.9789 
                + df_moles['K2O']*94.2 + df_moles['Cr2O3']*151.9904 + df_moles['NiO']*74.6928 + XFeO*71.844 + XFe2O3*159.69
                )
        
        #New Wt% values
        Fe2O3_wt = 100*(XFe2O3*159.69)/sum_moles
        FeO_wt = 100*(XFeO*71.844)/sum_moles



        if df_normalize == False:
            df_out = df[comps_noFe].copy()
            df_out['Fe2O3_calc'] = Fe2O3_wt
            df_out['FeO_calc'] = FeO_wt
            df_out = df_out[comps_final]
            df_out['Old_Total'] = pd.DataFrame.sum(df[comps], axis=1)
            df_out['New_Total'] = pd.DataFrame.sum(df[comps_noFe], axis=1) + Fe2O3_wt + FeO_wt
            df_out['Fe3Fetot'] = 2*XFe2O3/(2*XFe2O3+XFeO)
            df_out['XFe2O3_XFeO'] = XFe2O3_XFeO


        if df_normalize == True:
            df_out = pd.DataFrame(data={'SiO2': 100*df_moles['SiO2']*60.084/sum_moles,
                                    'TiO2': 100*df_moles['TiO2']*79.866/sum_moles,
                                    'Al2O3':100*df_moles['Al2O3']*101.96/sum_moles,
                                    'Fe2O3_calc': Fe2O3_wt,
                                    'FeO_calc': FeO_wt,
                                    'MnO': 100*df_moles['MnO']*70.9374/sum_moles,
                                    'MgO': 100*df_moles['MgO']*40.3044/sum_moles,
                                    'CaO': 100*df_moles['CaO']*56.0774/sum_moles,
                                    'Na2O': 100*df_moles['Na2O']*61.9789/sum_moles,
                                    'K2O': 100*df_moles['K2O']*94.2/sum_moles,
                                    'P2O5':  100*df_moles['P2O5']*141.943/sum_moles,
                                    'Cr2O3': 100*df_moles['Cr2O3']*151.9904/sum_moles,
                                    'NiO': 100*df_moles['NiO']*74.6928/sum_moles
                                    })

            df_out = df_out[comps_final]
            df_out['total'] = df_out.sum(axis=1)
            df_out['Fe3Fetot'] = 2*XFe2O3/(2*XFe2O3+XFeO)
            df_out['XFe2O3_XFeO'] = XFe2O3_XFeO





    elif model == "Kress_1991_lowp":

        XFe2O3_XFeO = np.exp(0.196*np.log(fO2)
                + 11492/T_K
                - 6.675
                - 2.243*df_moles['Al2O3']
                - 1.828*df_moles['FeO']
                + 3.201*df_moles['CaO']
                + 5.854*df_moles['Na2O']
                + 6.215*df_moles['K2O']
                )
        
        #Calculate new mole fractions
        XFe2O3 = XFe2O3_XFeO*df_moles.FeO/(2*XFe2O3_XFeO + 1)
        XFeO = df_moles['FeO'] - 2*XFe2O3


        sum_moles = (df_moles['SiO2']*60.084 + df_moles['TiO2']*79.866 + df_moles['Al2O3']*101.96 + df_moles['MnO']*70.9374 
                + df_moles['MgO']*40.3044 + df_moles['CaO']*56.0774 + df_moles['P2O5']*141.943 + df_moles['Na2O']*61.9789 
                + df_moles['K2O']*94.2 + df_moles['Cr2O3']*151.9904 + df_moles['NiO']*74.6928 + XFeO*71.844 + XFe2O3*159.69
                )
        
        #New Wt% values
        Fe2O3_wt = 100*(XFe2O3*159.69)/sum_moles
        FeO_wt = 100*(XFeO*71.844)/sum_moles



        if df_normalize == False:
            df_out = df[comps_noFe].copy()
            df_out['Fe2O3_calc'] = Fe2O3_wt
            df_out['FeO_calc'] = FeO_wt
            df_out = df_out[comps_final]
            df_out['Old_Total'] = pd.DataFrame.sum(df[comps], axis=1)
            df_out['New_Total'] = pd.DataFrame.sum(df[comps_noFe], axis=1) + Fe2O3_wt + FeO_wt
            df_out['Fe3Fetot'] = 2*XFe2O3/(2*XFe2O3+XFeO)
            df_out['XFe2O3_XFeO'] = XFe2O3_XFeO


        if df_normalize == True:
            df_out = pd.DataFrame(data={'SiO2': 100*df_moles['SiO2']*60.084/sum_moles,
                                    'TiO2': 100*df_moles['TiO2']*79.866/sum_moles,
                                    'Al2O3':100*df_moles['Al2O3']*101.96/sum_moles,
                                    'Fe2O3_calc': Fe2O3_wt,
                                    'FeO_calc': FeO_wt,
                                    'MnO': 100*df_moles['MnO']*70.9374/sum_moles,
                                    'MgO': 100*df_moles['MgO']*40.3044/sum_moles,
                                    'CaO': 100*df_moles['CaO']*56.0774/sum_moles,
                                    'Na2O': 100*df_moles['Na2O']*61.9789/sum_moles,
                                    'K2O': 100*df_moles['K2O']*94.2/sum_moles,
                                    'P2O5':  100*df_moles['P2O5']*141.943/sum_moles,
                                    'Cr2O3': 100*df_moles['Cr2O3']*151.9904/sum_moles,
                                    'NiO': 100*df_moles['NiO']*74.6928/sum_moles
                                    })

            df_out = df_out[comps_final]
            df_out['total'] = df_out.sum(axis=1)
            df_out['Fe3Fetot'] = 2*XFe2O3/(2*XFe2O3+XFeO)
            df_out['XFe2O3_XFeO'] = XFe2O3_XFeO
            




    elif model == "Jayasuriya_2004":
        XFe2O3_XFeO = np.exp(0.1967*np.log(fO2)
                                        + 12420/T_K
                                        - 7.054
                                        - 0.487*df_moles['MgO']
                                        + 2.201*df_moles['CaO']
                                        + 6.610*df_moles['Na2O']
                                        + 8.214*df_moles['K2O']
                                        - 3.781*df_moles['Al2O3']
                                        - 62.79*df_moles['P2O5']
                                        + 1.377*df_moles['FeO']
                                        )

        #Calculate new mole fractions
        XFe2O3 = XFe2O3_XFeO*df_moles.FeO/(2*XFe2O3_XFeO + 1)
        XFeO = df_moles['FeO'] - 2*XFe2O3


        sum_moles = (df_moles['SiO2']*60.084 + df_moles['TiO2']*79.866 + df_moles['Al2O3']*101.96 + df_moles['MnO']*70.9374 
                + df_moles['MgO']*40.3044 + df_moles['CaO']*56.0774 + df_moles['P2O5']*141.943 + df_moles['Na2O']*61.9789 
                + df_moles['K2O']*94.2 + df_moles['Cr2O3']*151.9904 + df_moles['NiO']*74.6928 + XFeO*71.844 + XFe2O3*159.69
                )
        
        #New Wt% values
        Fe2O3_wt = 100*(XFe2O3*159.69)/sum_moles
        FeO_wt = 100*(XFeO*71.844)/sum_moles



        if df_normalize == False:
            df_out = df[comps_noFe].copy()
            df_out['Fe2O3_calc'] = Fe2O3_wt
            df_out['FeO_calc'] = FeO_wt
            df_out = df_out[comps_final]
            df_out['Old_Total'] = pd.DataFrame.sum(df[comps], axis=1)
            df_out['New_Total'] = pd.DataFrame.sum(df[comps_noFe], axis=1) + Fe2O3_wt + FeO_wt
            df_out['Fe3Fetot'] = 2*XFe2O3/(2*XFe2O3+XFeO)
            df_out['XFe2O3_XFeO'] = XFe2O3_XFeO


        if df_normalize == True:
            df_out = pd.DataFrame(data={'SiO2': 100*df_moles['SiO2']*60.084/sum_moles,
                                    'TiO2': 100*df_moles['TiO2']*79.866/sum_moles,
                                    'Al2O3':100*df_moles['Al2O3']*101.96/sum_moles,
                                    'Fe2O3_calc': Fe2O3_wt,
                                    'FeO_calc': FeO_wt,
                                    'MnO': 100*df_moles['MnO']*70.9374/sum_moles,
                                    'MgO': 100*df_moles['MgO']*40.3044/sum_moles,
                                    'CaO': 100*df_moles['CaO']*56.0774/sum_moles,
                                    'Na2O': 100*df_moles['Na2O']*61.9789/sum_moles,
                                    'K2O': 100*df_moles['K2O']*94.2/sum_moles,
                                    'P2O5':  100*df_moles['P2O5']*141.943/sum_moles,
                                    'Cr2O3': 100*df_moles['Cr2O3']*151.9904/sum_moles,
                                    'NiO': 100*df_moles['NiO']*74.6928/sum_moles
                                    })

            df_out = df_out[comps_final]
            df_out['total'] = df_out.sum(axis=1)
            df_out['Fe3Fetot'] = 2*XFe2O3/(2*XFe2O3+XFeO)
            df_out['XFe2O3_XFeO'] = XFe2O3_XFeO

        



    elif model == "Putirka_2016":
        XFe2O3_XFeO = np.exp(-6.53 + 10813.8/T_K
        + 0.19*np.log(fO2)
        + 12.4*(df_moles['Na2O'] + df_moles['K2O'])
        - 3.44*(df_moles['Al2O3']/(df_moles['Al2O3']+df_moles['SiO2']))
        + 4.15*df_moles['CaO']
        )

        #Calculate new mole fractions
        XFe2O3 = XFe2O3_XFeO*df_moles.FeO/(2*XFe2O3_XFeO + 1)
        XFeO = df_moles['FeO'] - 2*XFe2O3


        sum_moles = (df_moles['SiO2']*60.084 + df_moles['TiO2']*79.866 + df_moles['Al2O3']*101.96 + df_moles['MnO']*70.9374 
                + df_moles['MgO']*40.3044 + df_moles['CaO']*56.0774 + df_moles['P2O5']*141.943 + df_moles['Na2O']*61.9789 
                + df_moles['K2O']*94.2 + df_moles['Cr2O3']*151.9904 + df_moles['NiO']*74.6928 + XFeO*71.844 + XFe2O3*159.69
                )
        
        #New Wt% values
        Fe2O3_wt = 100*(XFe2O3*159.69)/sum_moles
        FeO_wt = 100*(XFeO*71.844)/sum_moles



        if df_normalize == False:
            df_out = df[comps_noFe].copy()
            df_out['Fe2O3_calc'] = Fe2O3_wt
            df_out['FeO_calc'] = FeO_wt
            df_out = df_out[comps_final]
            df_out['Old_Total'] = pd.DataFrame.sum(df[comps], axis=1)
            df_out['New_Total'] = pd.DataFrame.sum(df[comps_noFe], axis=1) + Fe2O3_wt + FeO_wt
            df_out['Fe3Fetot'] = 2*XFe2O3/(2*XFe2O3+XFeO)
            df_out['XFe2O3_XFeO'] = XFe2O3_XFeO


        if df_normalize == True:
            df_out = pd.DataFrame(data={'SiO2': 100*df_moles['SiO2']*60.084/sum_moles,
                                    'TiO2': 100*df_moles['TiO2']*79.866/sum_moles,
                                    'Al2O3':100*df_moles['Al2O3']*101.96/sum_moles,
                                    'Fe2O3_calc': Fe2O3_wt,
                                    'FeO_calc': FeO_wt,
                                    'MnO': 100*df_moles['MnO']*70.9374/sum_moles,
                                    'MgO': 100*df_moles['MgO']*40.3044/sum_moles,
                                    'CaO': 100*df_moles['CaO']*56.0774/sum_moles,
                                    'Na2O': 100*df_moles['Na2O']*61.9789/sum_moles,
                                    'K2O': 100*df_moles['K2O']*94.2/sum_moles,
                                    'P2O5':  100*df_moles['P2O5']*141.943/sum_moles,
                                    'Cr2O3': 100*df_moles['Cr2O3']*151.9904/sum_moles,
                                    'NiO': 100*df_moles['NiO']*74.6928/sum_moles
                                    })

            df_out = df_out[comps_final]
            df_out['total'] = df_out.sum(axis=1)
            df_out['Fe3Fetot'] = 2*XFe2O3/(2*XFe2O3+XFeO)
            df_out['XFe2O3_XFeO'] = XFe2O3_XFeO
    




    elif model == "ONeill_2018":
        Fe3_Fe2 = 10**(0.25*(np.log10(fO2)-(8.58-(25050/T_K)))
                                  - 1.36
                                  + 2.4*df_moles_sc['CaO']
                                  + 2.0*df_moles_sc['Na2O']
                                  + 3.7*df_moles_sc['K2O']
                                  - 2.4*df_moles_sc['P2O5']
                                )
        Fe3_Fetot = 1/((1/Fe3_Fe2)+1)
        Fe2O3_FeO = Fe3_Fe2*(159.69/(2*71.844))

        Fe2O3_wt = (Fe2O3_FeO*df.FeO/(0.8998*Fe2O3_FeO + 1))
        FeO_wt = df.FeO - 0.8998*Fe2O3_wt

        if df_normalize == False:
            df_out = df[comps_noFe].copy()
            df_out['Fe2O3_calc'] = Fe2O3_wt
            df_out['FeO_calc'] = FeO_wt
            df_out = df_out[comps_final]
            df_out['Old_Total'] = pd.DataFrame.sum(df[comps], axis=1)
            df_out['New_Total'] = pd.DataFrame.sum(df[comps_noFe], axis=1) + Fe2O3_wt + FeO_wt
            df_out['Fe3Fetot'] = Fe3_Fetot
            df_out['Fe2O3_FeO'] = Fe2O3_FeO
            df_out['Fe3_Fe2'] = Fe3_Fe2


        if df_normalize == True:
            df_out['Fe2O3_calc'] = Fe2O3_wt
            df_out['FeO_calc'] = FeO_wt
            df_out[comps_noFe] = (normalize(df[comps_noFe].copy(), comps_noFe)).multiply(((100-df_out['Fe2O3_calc']-df_out['FeO_calc'])/100), axis=0)
            df_out = df_out[comps_final]
            df_out['total'] = df_out.sum(axis=1)
            df_out['Fe3Fetot'] = Fe3_Fetot
            df_out['Fe2O3_FeO'] = Fe2O3_FeO
            df_out['Fe3_Fe2'] = Fe3_Fe2





    elif model == "Borisov_2018":
        XFeO15_XFeO = 10**(0.207*np.log10(fO2)
                                 + 4633.3/T_K
                                 - 0.445*df_moles_borisov['SiO2']
                                 - 0.900*df_moles_borisov['TiO2']
                                 + 1.532*df_moles_borisov['MgO']
                                 + 0.314*df_moles_borisov['CaO']
                                 + 2.030*df_moles_borisov['Na2O']
                                 + 3.355*df_moles_borisov['K2O']
                                 - 4.851*df_moles_borisov['P2O5']
                                 -3.081*df_moles_borisov['SiO2']*df_moles_borisov['Al2O3']
                                 -4.370*df_moles_borisov['SiO2']*df_moles_borisov['MgO']
                                 -1.852
                                )
        XFeO15= (XFeO15_XFeO*df_moles_borisov.FeO/(XFeO15_XFeO + 1))
        XFeO = df_moles['FeO']/(1+XFeO15_XFeO)

        sum_moles = (df_moles_borisov['SiO2']*60.084 + df_moles_borisov['TiO2']*79.866 + df_moles_borisov['Al2O3']*101.96 
        + df_moles_borisov['MgO']*40.3044 + df_moles_borisov['CaO']*56.0774 + df_moles_borisov['P2O5']*141.943 + df_moles_borisov['Na2O']*61.9789 
        + df_moles_borisov['K2O']*94.2 + XFeO*71.844 + XFeO15*159.69*0.5
        )

        #New Wt% values
        Fe2O3_wt = 100*(XFeO15*159.69*0.5)/sum_moles
        FeO_wt = 100*(XFeO*71.844)/sum_moles


        if df_normalize == False:
            df_out = df[comps_noFe].copy()
            df_out['Fe2O3_calc'] = Fe2O3_wt
            df_out['FeO_calc'] = FeO_wt
            df_out = df_out[comps_final]
            df_out['Old_Total'] = pd.DataFrame.sum(df[comps], axis=1)
            df_out['New_Total'] = pd.DataFrame.sum(df[comps_noFe], axis=1) + Fe2O3_wt + FeO_wt
            df_out['Fe3Fetot'] = XFeO15/(XFeO15+XFeO)
            df_out['XFeO15_XFeO'] = XFeO15_XFeO


        if df_normalize == True:
            df_out = pd.DataFrame(data={'SiO2': 100*df_moles_borisov['SiO2']*60.084/sum_moles,
                                    'TiO2': 100*df_moles_borisov['TiO2']*79.866/sum_moles,
                                    'Al2O3':100*df_moles_borisov['Al2O3']*101.96/sum_moles,
                                    'Fe2O3': Fe2O3_wt,
                                    'FeO': FeO_wt,
                                    'MnO': df['MnO'],
                                    'MgO': 100*df_moles_borisov['MgO']*40.3044/sum_moles,
                                    'CaO': 100*df_moles_borisov['CaO']*56.0774/sum_moles,
                                    'Na2O': 100*df_moles_borisov['Na2O']*61.9789/sum_moles,
                                    'K2O': 100*df_moles_borisov['K2O']*94.2/sum_moles,
                                    'P2O5':  100*df_moles_borisov['P2O5']*141.943/sum_moles,
                                    'Cr2O3': df['Cr2O3'],
                                    'NiO': df['NiO']
                                    })

            df_out[comps_normalize] = (normalize(df_out[comps_normalize].copy(), comps_normalize))
            df_out.rename(columns={'FeO':'FeO_calc', 'Fe2O3':'Fe2O3_calc'}, inplace=True)
            df_out = df_out[comps_final]
            df_out['total'] = df_out.sum(axis=1)
            df_out['Fe3Fetot'] = XFeO15/(XFeO15+XFeO)
            df_out['XFeO15_XFeO'] = XFeO15_XFeO





    else:
        raise NotImplementedError(
            "Model {} not recognised.".format(model)
        )

    return df_out.join(index)





def monte_carlo(liq_comps, liq_comps_std, sample_row, mc):
    """
    Monte Carlo calculation model (Eason's code) for single sample
    Parameters:
    --------
    comp: :class: `pandas.Dataframe`
    std: :class: `pandas.Dataframe`
    sample_row: :class: `numpy.array`
        Row number of sample MC simulations are being performed for - will be used to rematch sample names after calculation
    mc: :class: `numpy.array`
        Number of monte carlo simulations to perform on selected sample
    Returns:
    ---------
    df: :class: `pandas.Dataframe`
        expanded dataframe after monte carlo simulation
    """
    liq_comps = liq_comps.copy()
    liq_comps_std = liq_comps_std.copy()
    index = sample_row
    index = np.ones(mc) * index                                             #Creates array with the row number of your selected sample, to allow names to be matched back later on

    
    liq_comps_mc =  liq_comps.iloc[sample_row,:]                            #Selects row of current sample and all columns from liq_comps dataframe
    liq_comps_std_mc = liq_comps_std.iloc[sample_row,:]                     #Selects row of current sample and all columns from liq_comps_std dataframe
    liq_comps_mc = np.ones((mc,1)) * liq_comps_mc.to_numpy()                #Create np array with mc rows and 1 column, then multiply by sample liquid composition to duplicate your sample composition and create a new array with mc rows
    liq_comps_std_mc = np.random.normal(0,1,(mc, len(liq_comps_std_mc))) * liq_comps_std_mc.to_numpy()                  #Generate a normally distributed array with a mean of 0 and std of 1, with mc rows and columns=number of oxides. Mulptiply by the std of your sample to create array of errors
    liq_comps_mc = liq_comps_mc + liq_comps_std_mc                          #Add your liquid composition to the error array to generate array of mc liquids that have compositions that vary within the stated error bounds
    df = pd.DataFrame(liq_comps_mc)                                         #Turn the array back into a dataframe
    df.index = index                                                        #Set the index for this new dataframe to be the same as the row number for the original sample (so all mc rows will have the same index number, to allow sample names to be matched back up)

    return df





def monte_carlo_iter(df_input, mc):
    """
    Monte Carlo calculation for dataframe of multiple samples
    Parameters:
    --------
    df_input: :class: `pandas.Dataframe`
        Dataframe containing sample index, temperature, pressure, logfO2, liquid compositions, standard deviations
    mc: :class: `numpy.array`
        Number of monte carlo simulations to perform on selected sample
    Returns:
    ---------
    df_out: :class: `pandas.Dataframe`
        expanded dataframe after monte carlo simulation for all samples
    """
    df_input = df_input.copy()
    mc = mc    

    columns = ['SiO2', 'TiO2', 'Al2O3', 'FeO', 'MnO', 'MgO',
           'CaO', 'Na2O', 'K2O', 'P2O5', 'NiO', 'Cr2O3']
    std_col = ['SiO2_std','TiO2_std', 'Al2O3_std', 'FeO_std', 'MnO_std', 'MgO_std', 
           'CaO_std','Na2O_std', 'K2O_std', 'P2O5_std', 'NiO_std', 'Cr2O3_std']

    liq_comps = df_input[columns].fillna(0)
    liq_comps_std = df_input[std_col].fillna(0)

    df_out = pd.DataFrame()
    for sample_row in range(len(liq_comps)):                                #Performs mc simulations on each sample in liq_comps df and combines into df_out
        df_iter = monte_carlo(liq_comps, liq_comps_std, sample_row, mc)
        df_out = pd.concat([df_out, df_iter])

    df_out.columns = liq_comps.columns                                          #Renames monte carlo columns to oxide headings
    df_out[df_out<0] = 0  	                                                    #Any negative oxide values will be set to 0
    df_out['sample_index'] = df_input.loc[df_input.index.repeat(mc)]['sample_index']  # match index name
    df_out['temperature'] = df_input.loc[df_input.index.repeat(mc)]['temperature']
    df_out['pressure'] = df_input.loc[df_input.index.repeat(mc)]['pressure']    #match exp pressure
    df_out['log_fO2'] = df_input.loc[df_input.index.repeat(mc)]['log_fO2']      #match exp temp
    df_out.reset_index(inplace=True)
                       
    return df_out