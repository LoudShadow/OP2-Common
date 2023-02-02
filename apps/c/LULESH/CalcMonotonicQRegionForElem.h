// Stuff needed for boundary conditions
// 2 BCs on each of 6 hexahedral faces (12 bits)

#define XI_M        0x00007
#define XI_M_SYMM   0x00001
#define XI_M_FREE   0x00002
#define XI_M_COMM   0x00004

// #define XI_P        0x00038
// #define XI_P_SYMM   0x00008
// #define XI_P_FREE   0x00010
// #define XI_P_COMM   0x00020

// #define ETA_M       0x001c0
// #define ETA_M_SYMM  0x00040
// #define ETA_M_FREE  0x00080
// #define ETA_M_COMM  0x00100

// #define ETA_P       0x00e00
// #define ETA_P_SYMM  0x00200
// #define ETA_P_FREE  0x00400
// #define ETA_P_COMM  0x00800

// #define ZETA_M      0x07000
// #define ZETA_M_SYMM 0x01000
// #define ZETA_M_FREE 0x02000
// #define ZETA_M_COMM 0x04000

// #define ZETA_P      0x38000
// #define ZETA_P_SYMM 0x08000
// #define ZETA_P_FREE 0x10000
// #define ZETA_P_COMM 0x20000

inline void CalcMonotonicQRegionForElem(
    const double *delv_xi, const double *delv_xi_lxim, const double *delv_xi_lxip,
    const double *delv_eta, const double *delv_eta_letam, const double *delv_eta_letap,
    const double *delv_zeta, const double *delv_zeta_lzetam, const double *delv_zeta_lzetap,
    const double *delx_xi, const double *delx_eta, const double *delx_zeta,
    const int *elemBC,
    const double *m_vdov,
    double *qq, double *ql,
    const double *elemMass, const double *volo, const double *vnew
){
    double qlin, qquad ;
    double phixi, phieta, phizeta ;
    double delvm = 0.0, delvp =0.0;
    int bcMask = elemBC[0] ;

    double norm = double(1.) / (delv_xi[0]+ m_ptiny ) ;

    switch (bcMask & XI_M) { //XI_M
    case 0x00004: /* needs comm data */ //XI_M_COMM
    case 0:         delvm = delv_xi_lxim[0]; break ;
    case 0x00001: delvm = delv_xi[0] ;       break ; //XI_M_SYMM
    case 0x00002: delvm = double(0.0) ;      break ; //XI_M_FREE
    default:          //fprintf(stderr, "Error in switch at %s line %d\n",__FILE__, __LINE__);
        delvm = 0; /* ERROR - but quiets the compiler */
        break;
    }
    switch (bcMask & 0x00038) { //XI_P
    case 0x00020: /* needs comm data */ //XI_P_COMM
    case 0:         delvp = delv_xi_lxip[0] ; break ;
    case 0x00008: delvp = delv_xi[0] ;       break ; //XI_P_FREE
    case 0x00010: delvp = double(0.0) ;      break ; //XI_P_FREE
    default:          //fprintf(stderr, "Error in switch at %s line %d\n",__FILE__, __LINE__);
        delvp = 0; /* ERROR - but quiets the compiler */
        break;
    }

    delvm = delvm * norm ;
    delvp = delvp * norm ;

    phixi = double(.5) * ( delvm + delvp ) ;

    delvm *= m_monoq_limiter_mult ;
    delvp *= m_monoq_limiter_mult ;

    if ( delvm < phixi ) phixi = delvm ;
    if ( delvp < phixi ) phixi = delvp ;
    if ( phixi < double(0.)) phixi = double(0.) ;
    if ( phixi > m_monoq_max_slope) phixi = m_monoq_max_slope;

    /*  phieta     */
    norm = double(1.) / ( delv_eta[0] + m_ptiny ) ;

    switch (bcMask & 0x001c0) { //ETA_M
        case 0x00100: /* needs comm data */ // ETA_M_COMM
        case 0:          delvm = delv_eta_letam[0] ; break ;
        case 0x00040: delvm = delv_eta[0] ;        break ; // ETA_M_SYMM
        case 0x00080: delvm = double(0.0) ;        break ; // ETA_M_FREE
        default:          //fprintf(stderr, "Error in switch at %s line %d\n",__FILE__, __LINE__);
        delvm = 0; /* ERROR - but quiets the compiler */
        break;
    }
    switch (bcMask & 0x00e00) { //ETA_P
        case 0x00800: /* needs comm data */  //ETA_P_COMM
        case 0:          delvp = delv_eta_letap[0] ; break ;
        case 0x00200: delvp = delv_eta[0] ;        break ; // ETA_P_SYMM
        case 0x00400: delvp = double(0.0) ;        break ; // ETA_P_FREE
        default:          //fprintf(stderr, "Error in switch at %s line %d\n",__FILE__, __LINE__);
        delvp = 0; /* ERROR - but quiets the compiler */
        break;
    }

    delvm = delvm * norm ;
    delvp = delvp * norm ;

    phieta = double(.5) * ( delvm + delvp ) ;

    delvm *= m_monoq_limiter_mult ;
    delvp *= m_monoq_limiter_mult ;

    if ( delvm  < phieta ) phieta = delvm ;
    if ( delvp  < phieta ) phieta = delvp ;
    if ( phieta < double(0.)) phieta = double(0.) ;
    if ( phieta > m_monoq_max_slope)  phieta = m_monoq_max_slope;

    /*  phizeta     */
    // norm = double(1.) / ( domain.delv_zeta(ielem) + ptiny ) ;
    norm = double(1.) / ( delv_zeta[0] + m_ptiny ) ;

    switch (bcMask & 0x07000) { //ZETA_M
        case 0x04000: /* needs comm data */ // ZETA_M_COMM
        case 0:           delvm = delv_zeta_lzetam[0] ; break ;
        case 0x01000: delvm = delv_zeta[0] ;         break ; // ZETA_M_SYMM
        case 0x02000: delvm = double(0.0) ;          break ; // ZETA_M_FREE
        default:          //fprintf(stderr, "Error in switch at %s line %d\n",__FILE__, __LINE__);
        delvm = 0; /* ERROR - but quiets the compiler */
        break;
    }
    switch (bcMask & 0x38000) { //ZETA_P
        case 0x20000: /* needs comm data */ // ZETA_P_COMM
        case 0:           delvp = delv_zeta_lzetap[0] ; break ;
        case 0x08000: delvp = delv_zeta[0] ;         break ; // ZETA_P_SYMM
        case 0x10000: delvp = double(0.0) ;          break ; // ZETA_P_FREE
        default:          //fprintf(stderr, "Error in switch at %s line %d\n",__FILE__, __LINE__);
        delvp = 0; /* ERROR - but quiets the compiler */
        break;
    }

    delvm = delvm * norm ;
    delvp = delvp * norm ;

    phizeta = double(.5) * ( delvm + delvp ) ;

    delvm *= m_monoq_limiter_mult ;
    delvp *= m_monoq_limiter_mult ;

    if ( delvm   < phizeta ) phizeta = delvm ;
    if ( delvp   < phizeta ) phizeta = delvp ;
    if ( phizeta < double(0.)) phizeta = double(0.);
    if ( phizeta > m_monoq_max_slope  ) phizeta = m_monoq_max_slope;

    /* Remove length scale */
    if ( m_vdov[0] > double(0.) )  {
        qlin  = double(0.) ;
        qquad = double(0.) ;
    }
    else {
        double delvxxi   = delv_xi[0]   * delx_xi[0]   ;
        double delvxeta  = delv_eta[0]  * delx_eta[0]  ;
        double delvxzeta = delv_zeta[0] * delx_zeta[0] ;

        if ( delvxxi   > double(0.) ) delvxxi   = double(0.) ;
        if ( delvxeta  > double(0.) ) delvxeta  = double(0.) ;
        if ( delvxzeta > double(0.) ) delvxzeta = double(0.) ;

        double rho = elemMass[0] / (volo[0] * vnew[0]) ;

        qlin = -m_qlc_monoq * rho *
        (  delvxxi   * (double(1.) - phixi) +
            delvxeta  * (double(1.) - phieta) +
            delvxzeta * (double(1.) - phizeta)  ) ;

        qquad = m_qqc_monoq * rho *
        (  delvxxi*delvxxi     * (double(1.) - phixi*phixi) +
            delvxeta*delvxeta   * (double(1.) - phieta*phieta) +
            delvxzeta*delvxzeta * (double(1.) - phizeta*phizeta)  ) ;
    }

    qq[0] = qquad ;
    ql[0] = qlin  ;

}