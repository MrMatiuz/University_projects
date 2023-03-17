#define IGD(k,l)  ( (k)*lt+(l) )

/*
#define max(a,b)  (((a) > (b)) ? (a) : (b))
#define min(a,b)  (((a) < (b)) ? (a) : (b))
*/

#define G_A_M(s)   ( (s)*1.6+(1.-(s))*1.4 )

    int riemi (int lom );
    void  ff_df( float p );
    double riem_c(double c_l,double c_r,double spd);

    void inp_param(void);
    void memo_on(void);
    void inp_reg(void);
    double step(void);
    void new_gdf(int k,int l);
    void deriv0(void);
    void deriv1(void);
    void deriv(int k);
    void t_deriv(int k, int l);
    double cnvt(int j, int k, int l);
    
    double x_flow0(void);
    void x_flow1(void);
    double x_flowk( int k );
    double  y_flow(int k,int l);
    double  y_flow1(int k);
    double  y_flow0(int k);
    void out_reg(void);
    double lim(double a,double b);
    
    
    
    
     
    
    
    
