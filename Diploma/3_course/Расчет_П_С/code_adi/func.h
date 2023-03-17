void inp_mesh(void);
void inp_gdf(void);
void derivs(void);
double step(void);
double flow_lr(void);
double flow_du(void);
void new_gdf(void);
double min_mod(double a, double b);
void out(void);
void out_tec(void);
void memo_on(void);

void riemi(double *par, int *alarm);


void fllrRusanov(int i, int j, double *frul, double *frur, double &db, double &ub, double &vb, double &pb, double *data);
void flduRusanov(int i, int j, double *frul, double *frur, double &db, double &ub, double &vb, double &pb, double *data);


/*double min(double,double)
double max(double,double)*/
