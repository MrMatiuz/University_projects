//ПРОГРАММА ДЛЯ НАХОЖДЕНИЯ ЛОКАЛЬНОГО МИНИМУМА ДЛЯ ФУНКЦИИ НЕСКОЛЬКИХ ПЕРЕМЕННЫХ ГРАДИЕНТНЫМ МЕТОДОМ И МЕТОДОМ НЬЮТОНА

#include <iostream>
#include <cmath> 
using namespace std;
const int n=3;
const double eps = 0.000001;
const double tay = 0.1*sqrt(eps);
const double betta = sqrt(eps);
//(pow(p[0],4)+pow(p[1],4)-pow(p[0]+p[1],2))
double func(double p[]){return (pow(p[0]-1,2)+0.1*pow(p[1]-2,2)+10*pow(p[2]-2,2));} // значение функции в точке

void derivative(double xyz[], double grad_x[]) // функция производной в точке xyz
{
	double norm=0;
	double res=0;
	for(int i=0; i<n; i++){
		xyz[i]+=tay; res=func(xyz);
		xyz[i]-=2*tay; res-=func(xyz);
		grad_x[i]=(res)/(2*tay);
		xyz[i]+=tay;
		//antigrad_x[i]=-grad_x[i];
		 //cout <<"производная" <<grad_x[i] << endl;
		norm+=pow(grad_x[i],2);
	}
	norm=sqrt(norm);
	for(int i=0; i<n; i++){
		grad_x[i]=grad_x[i]/norm;
		//cout<<"градиент"<<grad_x[i]<<endl;
	}
	//cout<<"----------"<<endl;
}

void second_derivative(double x[],double df[][n])// Функция нахождения матрицы вторых частных производных
{
    double res1=0;
    double res2=0;
	for(int i=0; i<n; i++)
	{
	    for(int j=0;j<n;j++)
	    {
	    	x[i]+=tay; res1=func(x);
	    	x[j]-=tay; res1-=func(x);
	    	x[i]-=tay; res1+=func(x);
	    	x[j]+=tay; res1-=func(x);
	    	res1=res1/pow(tay,2);

	    	x[j]+=tay; res2=func(x);
	    	x[i]-=tay; res2-=func(x);
	    	x[j]-=tay; res2+=func(x);
	    	x[i]+=tay; res2-=func(x);
	    	res2=res2/pow(tay,2);
	    	
	    	res2=(res1+res2)/2;
			df[i][j]=res2;
    		//cout << df[i][j] << " ";
	    }
	    //cout<<endl;
	}
}

void inverse_(double b[][n])// Функция подсчета обратной матрицы
{
   double a;
   double determinant = 0;
   for(int i=0;i<n;i++)
      determinant=determinant+(b[0][i]*(b[1][(i+1)%n]*b[2][(i+2)%n]-b[1][(i+2)%n]*b[2][(i+1)%n])); //finding determinant
   if(n==3){
      for(int i=0;i<n;i++)
      {
         for(int j=0;j<n;j++)
            b[i][j]=((b[(j+1)%n][(i+1)%n]*b[(j+2)%n][(i+2)%n])-(b[(j+1)%n][(i+2)%n]*b[(j+2)%n][(i+1)%n]))/ determinant;
      }
   }else{
		determinant=b[0][0]*b[1][1]-b[0][1]*b[1][0];
		a=b[0][0];
		b[0][0]=b[1][1]/determinant;
		b[0][1]=-b[0][1]/determinant;
		b[1][0]=-b[1][0]/determinant;
		b[1][1]=a/determinant;
    	}
}

void multi(double a[][n], double b[], double c[]) // Произведение матрицы NxN на Nx1; В вектор (с) записывается результат
{
	for(int i=0;i<n;i++)
	{
		c[i]=0;
		for(int j=0;j<n;j++)
			c[i]+=a[i][j]*b[j];
		//cout<<c[i]<<endl;
	}
}

void anti(double a[], double b[]) // функция смены знака
{
	for(int i=0; i<n; i++)
		a[i]=-b[i];
}

void show_arr(double a[]) // Вывод матрицы
{
	for(int i=0;i<n;i++)
		cout<<a[i]<<" ";
	cout<<endl;
}

void show_arr_2(double a[][n])
{
	for(int i=0; i<n; i++)
	{
	    for(int j=0;j<n;j++)
	    	cout<<a[i][j]<< " ";
	    cout<<endl;
	}
}

double step(double a[], double agr[], double l) // шаг градиентного метода
{
	double r[n];
	for(int i=0; i<n; i++)
		r[i]=a[i]+l*agr[i];
	return func(r);
}

double half(double a[], double agr[]) // Способ деления отрезка пополам для градиентного метода
{
	double check;
	int flag=0;
	double t=0;
	double a_n=0, b_n=0.01;
	double l1, l2;
	while( (b_n-a_n-eps)/pow(2,t+1 ) >= eps/2 )
	{
		//cout<<"t ";
		flag=0;
		t++;
		l1 = (a_n+b_n)/2 - eps/2; l2 = (a_n+b_n)/2 + eps/2;
		if(step(a,agr,l1) <= step(a,agr,l2))
			{b_n = l2; check=b_n;}
		else
		{
			a_n = l1;
			flag=1;
			check=a_n;
		}
	}
	if(flag==0) {/*cout<<"1  "<<a_n<<endl;*/ return a_n;}
	else {/*cout<<"2   "<<b_n<<endl;*/return b_n;}
}

double tangent(double a[], double agr[]) // Способ касательных для градиентного шага
{
	double a_n=0, b_n=0.01, c_n, flag;
	double l;
	int k=0;
	int m=0;
	double res=0;
	double ai=0, bi=0, ci=0;
	do{
		a_n+=tay; res=step(a,agr,a_n); //ВЫЧИСЛИЯЕМ ПРОИЗВОДНУЮ Ф-ИИ ФИ В ТОЧКЕ a_n
		a_n-=2*tay; res-=step(a,agr,a_n);
		ai=(res)/(2*tay);
		a_n+=tay;
		cout <<" a_n: "<<ai<<endl;

		b_n+=tay; res=step(a,agr,b_n); //ВЫЧИСЛИЯЕМ ПРОИЗВОДНУЮ Ф-ИИ ФИ В ТОЧКЕ b_n
		b_n-=2*tay; res-=step(a,agr,b_n);
		bi=(res)/(2*tay);
		b_n+=tay;
		cout <<" b_n: "<<bi<<endl;

		if(ai*bi>=0){k=1; break;}
		c_n=(step(a,agr,a_n)-step(a,agr,b_n)+bi*b_n-ai*a_n)/(bi-ai);

		c_n+=tay; res=step(a,agr,c_n); //ВЫЧИСЛИЯЕМ ПРОИЗВОДНУЮ Ф-ИИ ФИ В ТОЧКЕ c_n
		c_n-=2*tay; res-=step(a,agr,c_n);
		ci=(res)/(2*tay);
		c_n+=tay;

		if(flag!=c_n){
			if(ci<0) {a_n=c_n;}
			if(ci>0) {b_n=c_n;}
		}else {
			cout<<"exit";break;}
		
		flag=c_n;
		m++;
	//cout <<"точка с_n: "<<c_n<<endl;
	//cout <<"производная в точке с_n: "<<ci<<endl;
	//cout <<"--------------------------------------"<<endl;
	}while(abs(ci)>=eps || m==20);
	if(k==1){
		if(ai>=0 && bi>=0) return a_n;
		if(ai<=0 && bi<=0) return b_n;
	}
	return c_n;

}

double norm_subtraction(double a[], double b[]) // Норма разности двух векторов в евклидовом пространстве
{
	double res=0;
	double c[n];
	for(int i=0;i<n;i++)
	{
		c[i]=a[i]-b[i];
		res=res+pow(c[i],2);
	}
	return sqrt(res);
}

int if_(double a[],double b[],double coef) // Проверка условий окончания счета
{
	int c=2;
	double null_[n];
	double gr[n];
	for(int i=0;i<n;i++)
		null_[i]=0;
	derivative(a,gr); // Градиент в точке X k+1-ой 
	//cout<<norm_subtraction(a,b)<<endl;
	//cout<<norm_subtraction(gr,null_)<<endl;
	if(norm_subtraction(a,b)<=coef  && abs(func(a)-func(b))<=coef)//&& norm_subtraction(gr,null_)<=coef
		c=1;
	return c;
}

void grad_method() // Градиентный метод!!
{
	double xyz_o[n]={1.1,2.1,2.1}; //первая точка в окрестности лок. минимума
	int k=0;
	int m=0;
	double f_x_k_p;
	double lamda;
	double grad_x[n];
	double antigrad_x[n];
	double x_k[n]; //точка Х k-я
	double x_k_p[n];
	for(int i=0;i<n;i++)
		x_k[i]=xyz_o[i];
	do
	{
		if(k!=0){for(int i=0; i<n; i++) x_k[i]=x_k_p[i];}
		k++;
		derivative(x_k,grad_x); // Находим градиент
		anti(antigrad_x, grad_x); // Берем вектор, противополножный градиенту
		lamda=half(x_k,antigrad_x); // Определяем лямбду
		//cout<<lamda<<endl;
		f_x_k_p=step(x_k,antigrad_x,lamda); // Значение функции в точке X k+1-ой
		//cout<<"we"<<endl;
		for(int i=0; i<n; i++){
			x_k_p[i]=x_k[i]+lamda*antigrad_x[i]; // Делаем шаг градиентного метода
			//cout <<"Х-катая"<< x_k[i] << endl;
		}
		/*for(int i=0;i<n;i++)
		{
			if(abs(antigrad_x[i])<0.00000000000001)
				m++;
		}
		if(m!=0)
			break;
		*/
	}
	while(if_(x_k_p,x_k,betta)!=1);
	cout << "Количество шагов: " <<k<< endl;
	cout << "Точка минимума для вашей функции по градиентному методу: ";
	show_arr(x_k_p);
}

void newton_method() // Метод Ньютона!!
{
	double xyz_o[n]={1.05,2.05,2.05}; //первая точка в окрестности лок. минимума
	int k=0;
	int flag=0;
	double f_x_k_p;
	double lamda;
	double c[n]={1,1,1}; 
	double grad_x[n];
	double antigrad_x[n];
	double d2f[n][n];
	double x_k[n]; //точка Х k-я
	double x_k_p[n];
	for(int i=0;i<n;i++)
		x_k[i]=xyz_o[i];
	do
	{
		flag++;
		if(k!=0){for(int i=0; i<n; i++) x_k[i]=x_k_p[i];}
		k++;
		derivative(x_k,grad_x); // Находим градиент
		anti(antigrad_x, grad_x); // Берем вектор, противополножный градиенту
		second_derivative(x_k,d2f);
		inverse_(d2f); // Обращаем матрицу вторых производных
		multi(d2f,antigrad_x,c); // Нашли произведение антриградиента и матрицы вторых производных
		//show_arr(c);
		//cout<<endl;
		lamda=tangent(x_k,c); // Определяем лямбду
		//cout<<"вышел"<<endl;
		//cout<<"лямбда: "<<lamda<<endl;
		f_x_k_p=step(x_k,c,lamda); // Значение функции в точке X k+1-ой
		for(int i=0; i<n; i++){
			x_k_p[i]=x_k[i]+lamda*c[i]; // Делаем шаг градиентного метода
			//cout << x_k[i] <<" " << endl;
		}
	}
	while(if_(x_k_p,x_k,eps)!=1 || flag==200);
	cout << "Количество шагов: " <<k<< endl;
	cout << "Точка минимума для вашей функции по методу Ньютона: ";
	show_arr(x_k_p);
}

int main(){
	grad_method();
	newton_method();
	return 0;
}



