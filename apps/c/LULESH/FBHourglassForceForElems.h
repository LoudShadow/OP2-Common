inline void FBHourglassForceForElems(
    const double *p_xd0, const double *p_xd1, const double *p_xd2, const double *p_xd3, const double *p_xd4, const double *p_xd5, const double *p_xd6, const double *p_xd7,
    const double *p_yd0, const double *p_yd1, const double *p_yd2, const double *p_yd3, const double *p_yd4, const double *p_yd5, const double *p_yd6, const double *p_yd7,
    const double *p_zd0, const double *p_zd1, const double *p_zd2, const double *p_zd3, const double *p_zd4, const double *p_zd5, const double *p_zd6, const double *p_zd7,
    double *p_fx0, double *p_fx1, double *p_fx2, double *p_fx3, double *p_fx4, double *p_fx5, double *p_fx6, double *p_fx7,
    double *p_fy0, double *p_fy1, double *p_fy2, double *p_fy3, double *p_fy4, double *p_fy5, double *p_fy6, double *p_fy7,
    double *p_fz0, double *p_fz1, double *p_fz2, double *p_fz3, double *p_fz4, double *p_fz5, double *p_fz6, double *p_fz7,
    const double *p_dvdx,
    const double *p_dvdy,
    const double *p_dvdz,
    const double *p_x8n,
    const double *p_y8n,
    const double *p_z8n,
    const double *p_determ,
    const double *p_ss,
    const double *p_elemMass
){

    double hgfx[8], hgfy[8], hgfz[8] ;

    double coefficient;

    double hourgam[8][4];

    double volinv=double(1.0)/p_determ[0];
    double volume13 ;


        for(int i1=0;i1<4;++i1){
        double hourmodx =
            p_x8n[0] * m_gamma_t[(i1*8)+0] + p_x8n[0+1] * m_gamma_t[(i1*8)+1] +
            p_x8n[0+2] * m_gamma_t[(i1*8)+2] + p_x8n[0+3] * m_gamma_t[(i1*8)+3] +
            p_x8n[0+4] * m_gamma_t[(i1*8)+4] + p_x8n[0+5] * m_gamma_t[(i1*8)+5] +
            p_x8n[0+6] * m_gamma_t[(i1*8)+6] + p_x8n[0+7] * m_gamma_t[(i1*8)+7];

        double hourmody =
            p_y8n[0] * m_gamma_t[(i1*8)+0] + p_y8n[0+1] * m_gamma_t[(i1*8)+1] +
            p_y8n[0+2] * m_gamma_t[(i1*8)+2] + p_y8n[0+3] * m_gamma_t[(i1*8)+3] +
            p_y8n[0+4] * m_gamma_t[(i1*8)+4] + p_y8n[0+5] * m_gamma_t[(i1*8)+5] +
            p_y8n[0+6] * m_gamma_t[(i1*8)+6] + p_y8n[0+7] * m_gamma_t[(i1*8)+7];

        double hourmodz =
            p_z8n[0] * m_gamma_t[(i1*8)+0] + p_z8n[0+1] * m_gamma_t[(i1*8)+1] +
            p_z8n[0+2] * m_gamma_t[(i1*8)+2] + p_z8n[0+3] * m_gamma_t[(i1*8)+3] +
            p_z8n[0+4] * m_gamma_t[(i1*8)+4] + p_z8n[0+5] * m_gamma_t[(i1*8)+5] +
            p_z8n[0+6] * m_gamma_t[(i1*8)+6] + p_z8n[0+7] * m_gamma_t[(i1*8)+7];

        hourgam[0][i1] = m_gamma_t[(i1*8)+0] -  volinv*(p_dvdx[0  ] * hourmodx +
                                                p_dvdy[0  ] * hourmody +
                                                p_dvdz[0  ] * hourmodz );

        hourgam[1][i1] = m_gamma_t[(i1*8)+1] -  volinv*(p_dvdx[0+1] * hourmodx +
                                                p_dvdy[0+1] * hourmody +
                                                p_dvdz[0+1] * hourmodz );

        hourgam[2][i1] = m_gamma_t[(i1*8)+2] -  volinv*(p_dvdx[0+2] * hourmodx +
                                                p_dvdy[0+2] * hourmody +
                                                p_dvdz[0+2] * hourmodz );

        hourgam[3][i1] = m_gamma_t[(i1*8)+3] -  volinv*(p_dvdx[0+3] * hourmodx +
                                                p_dvdy[0+3] * hourmody +
                                                p_dvdz[0+3] * hourmodz );

        hourgam[4][i1] = m_gamma_t[(i1*8)+4] -  volinv*(p_dvdx[0+4] * hourmodx +
                                                p_dvdy[0+4] * hourmody +
                                                p_dvdz[0+4] * hourmodz );

        hourgam[5][i1] = m_gamma_t[(i1*8)+5] -  volinv*(p_dvdx[0+5] * hourmodx +
                                                p_dvdy[0+5] * hourmody +
                                                p_dvdz[0+5] * hourmodz );

        hourgam[6][i1] = m_gamma_t[(i1*8)+6] -  volinv*(p_dvdx[0+6] * hourmodx +
                                                p_dvdy[0+6] * hourmody +
                                                p_dvdz[0+6] * hourmodz );

        hourgam[7][i1] = m_gamma_t[(i1*8)+7] -  volinv*(p_dvdx[0+7] * hourmodx +
                                                p_dvdy[0+7] * hourmody +
                                                p_dvdz[0+7] * hourmodz );
    }

    /* compute forces */
    /* store forces into h arrays (force arrays) */

    volume13=cbrt(p_determ[0]);

    coefficient = - m_hgcoef * double(0.01) * p_ss[0] * p_elemMass[0] / volume13;

    double hxx[4];
    for(int i = 0; i < 4; i++) {
        hxx[i] = hourgam[0][i] * p_xd0[0] + hourgam[1][i] * p_xd1[0] +
                hourgam[2][i] * p_xd2[0] + hourgam[3][i] * p_xd3[0] +
                hourgam[4][i] * p_xd4[0] + hourgam[5][i] * p_xd5[0] +
                hourgam[6][i] * p_xd6[0] + hourgam[7][i] * p_xd7[0];
    }
    for(int i = 0; i < 8; i++) {
        hgfx[i] = coefficient *
                    (hourgam[i][0] * hxx[0] + hourgam[i][1] * hxx[1] +
                    hourgam[i][2] * hxx[2] + hourgam[i][3] * hxx[3]);
    }
    for(int i = 0; i < 4; i++) {
        hxx[i] = hourgam[0][i] * p_yd0[0] + hourgam[1][i] * p_yd1[0] +
                hourgam[2][i] * p_yd2[0] + hourgam[3][i] * p_yd3[0] +
                hourgam[4][i] * p_yd4[0] + hourgam[5][i] * p_yd5[0] +
                hourgam[6][i] * p_yd6[0] + hourgam[7][i] * p_yd7[0];
    }
    for(int i = 0; i < 8; i++) {
        hgfy[i] = coefficient *
                    (hourgam[i][0] * hxx[0] + hourgam[i][1] * hxx[1] +
                    hourgam[i][2] * hxx[2] + hourgam[i][3] * hxx[3]);
    }
    for(int i = 0; i < 4; i++) {
        hxx[i] = hourgam[0][i] * p_zd0[0] + hourgam[1][i] * p_zd1[0] +
                hourgam[2][i] * p_zd2[0] + hourgam[3][i] * p_zd3[0] +
                hourgam[4][i] * p_zd4[0] + hourgam[5][i] * p_zd5[0] +
                hourgam[6][i] * p_zd6[0] + hourgam[7][i] * p_zd7[0];
    }
    for(int i = 0; i < 8; i++) {
        hgfz[i] = coefficient *
                    (hourgam[i][0] * hxx[0] + hourgam[i][1] * hxx[1] +
                    hourgam[i][2] * hxx[2] + hourgam[i][3] * hxx[3]);
    }

    p_fx0[0] += hgfx[0];
    p_fy0[0] += hgfy[0];
    p_fz0[0] += hgfz[0];

    p_fx1[0] += hgfx[1];
    p_fy1[0] += hgfy[1];
    p_fz1[0] += hgfz[1];

    p_fx2[0] += hgfx[2];
    p_fy2[0] += hgfy[2];
    p_fz2[0] += hgfz[2];

    p_fx3[0] += hgfx[3];
    p_fy3[0] += hgfy[3];
    p_fz3[0] += hgfz[3];

    p_fx4[0] += hgfx[4];
    p_fy4[0] += hgfy[4];
    p_fz4[0] += hgfz[4];

    p_fx5[0] += hgfx[5];
    p_fy5[0] += hgfy[5];
    p_fz5[0] += hgfz[5];

    p_fx6[0] += hgfx[6];
    p_fy6[0] += hgfy[6];
    p_fz6[0] += hgfz[6];

    p_fx7[0] += hgfx[7];
    p_fy7[0] += hgfy[7];
    p_fz7[0] += hgfz[7];
}