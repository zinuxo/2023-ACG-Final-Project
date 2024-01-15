#include<bits/stdc++.h>
using namespace std;

const int N=2e5+5,N2=105;

const float pi=acos(-1);
const float eps=1e-4;
const float g=10;
const float mu=0.003;
const float rho0=1e3;
const float k=40;
const float dt=1e-3;
// const float dt=4e-4;
const float h=0.04;
const float omega=0.05;
const float gamma=0.7;


inline float sqr(float x){return x*x;}
inline float cubic(float x){return x*x*x;}
struct vec3{
	float x,y,z;
	inline float&operator[](int i){return i==0?x:i==1?y:z;}
	inline float operator[](int i)const{return i==0?x:i==1?y:z;}
};
inline float dot(vec3 a,vec3 b){return a.x*b.x+a.y*b.y+a.z*b.z;}
inline vec3 operator-(vec3 a,vec3 b){return vec3{a.x-b.x,a.y-b.y,a.z-b.z};}
inline vec3 operator+(vec3 a,vec3 b){return vec3{a.x+b.x,a.y+b.y,a.z+b.z};}
inline vec3 operator*(vec3 a,float b){return vec3{a.x*b,a.y*b,a.z*b};}
inline float len2(vec3 a){return dot(a,a);}
inline float len(vec3 a){return sqrt(len2(a));}
inline vec3 norm(vec3 a){return a*(1.0/len(a));}
inline vec3 cross(vec3 a,vec3 b){return vec3{a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
inline float getW(vec3 r){
	float q=len(r)/h;
	if(q<=1){
		return 3/(2*pi*h*h*h)*(2.0/3-q*q+0.5*q*q*q);
	}else if(q<=2){
		return 3/(2*pi*h*h*h)*1.0/6*(2-q)*(2-q)*(2-q);
	}else return 0;
}
inline vec3 getWG(vec3 r){
	if(len(r)<eps)return vec3{0,0,0};
	float q=len(r)/h;
	vec3 res=r*(1.0/(len(r)*h));
	if(q<=1){
		return res*(3/(2*pi*h*h*h)*(-2*q+1.5*q*q));
	}else if(q<=2){
		return res*(3/(2*pi*h*h*h)*(-1.0/2*(2-q)*(2-q)));
	}else return vec3{0,0,0};
}

const vec3 BLUE={0.1, 0.6, 0.9};
const vec3 YELLOW={0.9, 0.5, 0.1};
const vec3 WHITE={1,1,1};
const vec3 GREEN={0.1,0.9,0.1};
const vec3 colors[4]={BLUE,YELLOW,WHITE,GREEN};
const vec3 lb=vec3{-h,-h,-h};
const vec3 rb=vec3{3+h,2+h,2+h};
const int M1=ceil((rb.x-lb.x)/h/2),M2=ceil((rb.y-lb.y)/h/2),M3=ceil((rb.z-lb.z)/h/2);

struct mat3{
	float a[3][3];
};
inline mat3 operator*(mat3 a,mat3 b){
	mat3 c;
	for(int i=0;i<3;++i)for(int j=0;j<3;++j){
		c.a[i][j]=0;
		for(int k=0;k<3;++k){
			c.a[i][j]+=a.a[i][k]*b.a[k][j];
		}
	}
	return c;
}
inline vec3 operator*(mat3 a,vec3 b){
	vec3 c;
	for(int i=0;i<3;++i){
		c[i]=0;
		for(int j=0;j<3;++j){
			c[i]+=a.a[i][j]*b[j];
		}
	}
	return c;
}
inline mat3 operator+(mat3 a,mat3 b){
	mat3 c;
	for(int i=0;i<3;++i)for(int j=0;j<3;++j){
		c.a[i][j]=a.a[i][j]+b.a[i][j];
	}
	return c;
}
inline mat3 operator*(mat3 a,float b){
	mat3 c;
	for(int i=0;i<3;++i)for(int j=0;j<3;++j){
		c.a[i][j]=a.a[i][j]*b;
	}
	return c;
}
inline mat3 operator-(mat3 a,mat3 b){
	mat3 c;
	for(int i=0;i<3;++i)for(int j=0;j<3;++j){
		c.a[i][j]=a.a[i][j]-b.a[i][j];
	}
	return c;
}
inline mat3 getinv(mat3 a){
	mat3 b;
	float det=a.a[0][0]*a.a[1][1]*a.a[2][2]+a.a[0][1]*a.a[1][2]*a.a[2][0]+a.a[0][2]*a.a[1][0]*a.a[2][1]
		-a.a[0][0]*a.a[1][2]*a.a[2][1]-a.a[0][1]*a.a[1][0]*a.a[2][2]-a.a[0][2]*a.a[1][1]*a.a[2][0];
	assert(det!=0);
	b.a[0][0]=(a.a[1][1]*a.a[2][2]-a.a[1][2]*a.a[2][1])/det;
	b.a[0][1]=(a.a[0][2]*a.a[2][1]-a.a[0][1]*a.a[2][2])/det;
	b.a[0][2]=(a.a[0][1]*a.a[1][2]-a.a[0][2]*a.a[1][1])/det;
	b.a[1][0]=(a.a[1][2]*a.a[2][0]-a.a[1][0]*a.a[2][2])/det;
	b.a[1][1]=(a.a[0][0]*a.a[2][2]-a.a[0][2]*a.a[2][0])/det;
	b.a[1][2]=(a.a[0][2]*a.a[1][0]-a.a[0][0]*a.a[1][2])/det;
	b.a[2][0]=(a.a[1][0]*a.a[2][1]-a.a[1][1]*a.a[2][0])/det;
	b.a[2][1]=(a.a[0][1]*a.a[2][0]-a.a[0][0]*a.a[2][1])/det;
	b.a[2][2]=(a.a[0][0]*a.a[1][1]-a.a[0][1]*a.a[1][0])/det;
	return b;
}
inline mat3 a2m(vec3 a){
	mat3 b;
	b.a[0][0]=0;
	b.a[0][1]=-a.z;
	b.a[0][2]=a.y;
	b.a[1][0]=a.z;
	b.a[1][1]=0;
	b.a[1][2]=-a.x;
	b.a[2][0]=-a.y;
	b.a[2][1]=a.x;
	b.a[2][2]=0;
	return b;
}
const mat3 I3=mat3{1,0,0,0,1,0,0,0,1};


struct Rigid{
	float mass;
	vec3 c,v,w;
	mat3 I,Iinv,trans;
	void output(){
		printf("%.8f\n",mass);
		printf("%.8f %.8f %.8f\n",c.x,c.y,c.z);
		for(int i=0;i<3;++i){
			for(int j=0;j<3;++j){
				printf("%.8f ",I.a[i][j]);
			}
			printf("\n");
		}
		for(int i=0;i<3;++i){
			for(int j=0;j<3;++j){
				printf("%.8f ",Iinv.a[i][j]);
			}
			printf("\n");
		}
	}
	void input(){
		scanf("%f",&mass);
		scanf("%f%f%f",&c.x,&c.y,&c.z);
		for(int i=0;i<3;++i){
			for(int j=0;j<3;++j){
				scanf("%f",&I.a[i][j]);
			}
		}
		for(int i=0;i<3;++i){
			for(int j=0;j<3;++j){
				scanf("%f",&Iinv.a[i][j]);
			}
		}
		trans=I3;
	}
}ri[N2];
int rcnt;

/*
	0 fluid
	1 static rigid
	2 dynamic rigid
*/
struct Particles{
	vec3 x[N],v[N];
	float mass[N],rho[N];
	int info[N];
	int n;
	void copy(Particles&rhs){
		memcpy(x+1,rhs.x+1,n*sizeof(vec3));
		memcpy(v+1,rhs.v+1,n*sizeof(vec3));
		memcpy(mass+1,rhs.mass+1,n*sizeof(float));
		memcpy(rho+1,rhs.rho+1,n*sizeof(float));
		memcpy(info+1,rhs.info+1,n*sizeof(int));
	}
}A,B;
inline int get_id(vec3 p){
	int x=floor((p.x-lb.x)/h/2);
	int y=floor((p.y-lb.y)/h/2);
	int z=floor((p.z-lb.z)/h/2);
	if(!(0<=x && x<M1))cerr<<p.x<<' '<<p.y<<' '<<p.z<<endl;
	assert(0<=x && x<M1);
	assert(0<=y && y<M2);
	assert(0<=z && z<M3);
	return x*M2*M3+y*M3+z;
}
inline mat3 getK(int i,int j){
	int R=A.info[i]>>8;
	return I3*(1.0/ri[R].mass)-a2m(A.x[i]-ri[R].c)*ri[R].Iinv*a2m(A.x[j]-ri[R].c);
}

int be[M1*M2*M3];
inline void radix_sort(){	
	int i;
	memset(be,0,sizeof be);
	for(i=1;i<=A.n;++i){
		be[get_id(A.x[i])]++;
	}
	int max=0;
	for(i=1;i<M1*M2*M3;++i){
		max=max>be[i]?max:be[i];
		be[i]+=be[i-1];
	}
	// cerr<<max<<endl;
	B.n=A.n;
	for(i=A.n;i>=1;--i){
		int id=be[get_id(A.x[i])]--;
		B.x[id]=A.x[i];
		B.v[id]=A.v[i];
		B.mass[id]=A.mass[i];
		B.rho[id]=A.rho[i];
		B.info[id]=A.info[i];		
	}
	A.copy(B);
}

namespace Iterating_State{
	float aii[N],s[N],p[N],lp[N];
	vec3 v2[N],fnonp[N],tmp_ap[N];
	vec3 FR_nonp[N2],TR_nonp[N2],FR[N2],TR[N2],vs[N],vr[N];
	bool in_contact[N];
	int tot0,tot1;
	float sum_err0,sum_err1;
}
inline void update_rigid(int i,vec3 f){
	using namespace Iterating_State;
	int id=A.info[i]>>8;
	FR[id]=FR[id]+f;
	TR[id]=TR[id]+cross(A.x[i]-ri[id].c,f);
}
template<typename T>inline void for_neighbor(){
	T local_updater;
	for(int i=1;i<=A.n;++i)if(local_updater.check(i)){
		local_updater.init(i);
		int x=floor((A.x[i].x-lb.x)/h/2),y=floor((A.x[i].y-lb.y)/h/2),z=floor((A.x[i].z-lb.z)/h/2);
		for(int dx=-1;dx<=1;++dx)for(int dy=-1;dy<=1;++dy)for(int dz=-1;dz<=1;++dz){
			int xx=x+dx,yy=y+dy,zz=z+dz;
			if(xx<0 || xx>=M1 || yy<0 || yy>=M2 || zz<0 || zz>=M3)continue;
			int id=xx*M2*M3+yy*M3+zz,start=be[id]+1,end=id==M1*M2*M3-1?A.n:be[id+1];
			for(int j=start;j<=end;++j)if(len(A.x[i]-A.x[j])<=2*h){
				local_updater.update(i,j);
			}
		}
		local_updater.finish(i);
	}
}
struct update_rho{
	float new_rho;
	inline bool check(int i){
		return (A.info[i]&15)!=1;
	}
	inline void init(int i){
		new_rho=0;
	}
	inline void update(int i,int j){
		using namespace Iterating_State;
		// if((A.info[i]&15)!=2 || ((A.info[j]&15)!=0))
		new_rho+=A.mass[j]*getW(A.x[i]-A.x[j]);
		if((A.info[i]&15)!=0 && (A.info[j]&15)!=0 && A.info[i]!=A.info[j])
			in_contact[i]=true;
	}
	inline void finish(int i){
		new_rho=max(new_rho,rho0);
		A.rho[i]=new_rho;
	}
}U1;
static int stepp;
struct update2{
	vec3 sum1,nonp;
	float sum2;
	inline bool check(int i){
		return (A.info[i]&15)==0;
	}
	inline void init(int i){
		sum1=vec3{0,0,0};
		sum2=0;
		nonp=vec3{0,0,-g*A.mass[i]};
	}
	inline void update(int i,int j){
		if(len(A.x[i]-A.x[j])<eps)return;
		sum1=sum1+getWG(A.x[i]-A.x[j])*A.mass[j];
		sum2=sum2+len2(getWG(A.x[i]-A.x[j]))*A.mass[j];
		vec3 f=(A.v[i]-A.v[j])*(mu*(A.mass[j]/A.rho[j])*2*
			len(getWG(A.x[i]-A.x[j]))/len(A.x[i]-A.x[j]));
		nonp=nonp-f;
		if((A.info[j]&15)==2)
			update_rigid(j,f);
	}
	inline void finish(int i){
		using namespace Iterating_State;
		aii[i]=-(len2(sum1)+sum2*A.mass[i])*sqr(dt)/sqr(A.rho[i]);
		fnonp[i]=nonp;
		v2[i]=A.v[i]+fnonp[i]*(dt/A.mass[i]);
	}
};
struct update3{
	float sum;
	inline bool check(int i){
		return (A.info[i]&15)==0;
	}
	inline void init(int i){
		sum=0;
	}
	inline void update(int i,int j){
		using namespace Iterating_State;
		sum=sum+A.mass[j]*dot(v2[i]-v2[j],getWG(A.x[i]-A.x[j]));
	}
	inline void finish(int i){
		using namespace Iterating_State;
		s[i]=rho0-A.rho[i]-sum*dt;
	}
};
struct update_fluid_pressure_force{
	vec3 sum;
	inline bool check(int i){
		return (A.info[i]&15)==0;
	}
	inline void init(int i){
		sum=vec3{0,0,0};
	}
	inline void update(int i,int j){
		using namespace Iterating_State;
		vec3 f=getWG(A.x[i]-A.x[j])*(A.mass[j]*(lp[i]/sqr(A.rho[i])+
			((A.info[j]&15)==0?lp[j]/sqr(A.rho[j]):0)));
		sum=sum-f;
		if((A.info[j]&15)==2)
			update_rigid(j,f);
	}
	inline void finish(int i){
		using namespace Iterating_State;
		tmp_ap[i]=sum;
	}
};
struct update_s{
	float sum;
	inline bool check(int i){
		using namespace Iterating_State;
		return in_contact[i];
	}
	inline void init(int i){
		sum=0;
	}
	inline void update(int i,int j){
		using namespace Iterating_State;
		// if(i!=j && A.info[i]==A.info[j]){
		if(i!=j && (A.info[j]&15)!=1){
			// sum=sum+dot(vs[j]-vs[i],getWG(A.x[i]-A.x[j]))*A.mass[j];
			sum=sum+dot(vs[j]-vs[i],getWG(A.x[i]-A.x[j]))*A.mass[j];
		}
	}
	inline void finish(int i){
		using namespace Iterating_State;
		s[i]=(rho0-A.rho[i])/dt+sum;		
	}
};
struct update_rigid_pressure_force{
	vec3 sum;
	inline bool check(int i){
		using namespace Iterating_State;
		return in_contact[i];
	}
	inline void init(int i){
		sum=vec3{0,0,0};
	}
	inline void update(int i,int j){
		using namespace Iterating_State;
		if(A.info[i]!=A.info[j] && (A.info[j]&15)!=0){
			vec3 f=getWG(A.x[i]-A.x[j])*(A.mass[j]*(lp[i]/sqr(A.rho[i])+
				((A.info[j]&15)==2?lp[j]/sqr(A.rho[j]):0)));
			sum=sum-f;
if(isnan(f.x))
	++i,--i;
		}
	}
	inline void finish(int i){
		using namespace Iterating_State;
		tmp_ap[i]=sum;
	}
};
bool flag;
struct update4{
	float sum=0;
	inline bool check(int i){
		using namespace Iterating_State;
		return (A.info[i]&15)==0 || in_contact[i];
	}
	inline void init(int i){
		sum=0;
	}
	inline void update(int i,int j){
		using namespace Iterating_State;
		if(in_contact[i]){
			if((A.info[j]&15)==0)return;
			// if(A.info[i]!=A.info[j])return;
			sum=sum-dot((A.info[i]==A.info[j]?vr[j]:vec3{0,0,0})-vr[i],getWG(A.x[i]-A.x[j]))*A.mass[j];
			++i,--i;
		}else{
			// if((A.info[j]&15)!=0)return;
			sum=sum+A.mass[j]*dot(tmp_ap[i]-((A.info[j]&15)==0?tmp_ap[j]:vec3{0,0,0}),getWG(A.x[i]-A.x[j]));
		}
	}
	inline void finish(int i){
		using namespace Iterating_State;
		if((A.info[i]&15)==0){
			sum*=sqr(dt);
		}
		if(in_contact[i] && s[i]==0)return;
		if((A.info[i]&15)==0 && flag){
			return;
		}
if(in_contact[i])
++i,--i;
if((s[i]-sum)!=0)
++i,--i;
		if(fabs(aii[i])>(in_contact[i]?5e-3:1e-5))
				p[i]=max(lp[i]+(in_contact[i]?0.1f:omega)/aii[i]*(s[i]-sum),0.0f);
			else p[i]=0;
		float err=fabs(s[i]-sum)/rho0;
		float&sum_err=(A.info[i]&15)==0?sum_err0:sum_err1;
		int&tot=(A.info[i]&15)==0?tot0:tot1;
		sum_err+=err;
		++tot;
	}
};

void step(){
	using namespace Iterating_State;
	radix_sort();
	memset(in_contact+1,0,A.n);
	for_neighbor<update_rho>();
	for(int i=1;i<=A.n;++i)v2[i]=tmp_ap[i]=vr[i]=vec3{0,0,0};
	for(int i=1;i<=A.n;++i)lp[i]=p[i]=aii[i]=0;
	for(int i=1;i<=rcnt;++i)FR[i]=vec3{0,0,-g*ri[i].mass},TR[i]=vec3{0,0,0};
cerr<<++stepp<<endl;
if(stepp==65)
++stepp,--stepp;
	for_neighbor<update2>();
	for_neighbor<update3>();
	memcpy(FR_nonp+1,FR+1,rcnt*sizeof(vec3));
	memcpy(TR_nonp+1,TR+1,rcnt*sizeof(vec3));
	sum_err0=0;
	tot0=0;
	sum_err1=0;
	tot1=0;
	flag=0;
	for(int i=1;i<=A.n;++i)if(in_contact[i]){
		static int id[999];int xb=0;
		int x=floor((A.x[i].x-lb.x)/h/2),y=floor((A.x[i].y-lb.y)/h/2),z=floor((A.x[i].z-lb.z)/h/2);
		for(int dx=-1;dx<=1;++dx)for(int dy=-1;dy<=1;++dy)for(int dz=-1;dz<=1;++dz){
			int xx=x+dx,yy=y+dy,zz=z+dz;
			if(xx<0 || xx>=M1 || yy<0 || yy>=M2 || zz<0 || zz>=M3)continue;
			int idd=xx*M2*M3+yy*M3+zz,start=be[idd]+1,end=idd==M1*M2*M3-1?A.n:be[idd+1];
			for(int j=start;j<=end;++j)if((A.info[j]&15)!=0 && len(A.x[i]-A.x[j])<=2*h){
				id[++xb]=j;
			}
		}
		vec3 sum_f=vec3{0,0,0};
		for(int jj=1;jj<=xb;++jj)if(A.info[i]!=A.info[id[jj]]){
			int j=id[jj];
			sum_f=sum_f-getWG(A.x[i]-A.x[j])*(A.mass[j]*(1/sqr(A.rho[i]))*A.mass[i]*dt);
		}
		float sum=0;
		for(int kk=1;kk<=xb;++kk){
			int k=id[kk];
			if(i==k)continue;
			// if(i==k || A.info[i]!=A.info[k])continue;
			sum=sum-dot((A.info[k]==A.info[i]?getK(k,i)*sum_f:vec3{0,0,0})
				-getK(i,i)*sum_f,getWG(A.x[i]-A.x[k]))*A.mass[k];
		}
		aii[i]=sum;
	}
	// for(int i=1;i<=A.n;++i)if(in_contact[i])cerr<<i<<' ';cerr<<endl;
	float last;
	for(int ite=0;ite<10;++ite){
		memcpy(FR+1,FR_nonp+1,rcnt*sizeof(vec3));
		memcpy(TR+1,TR_nonp+1,rcnt*sizeof(vec3));
		if(!flag)for_neighbor<update_fluid_pressure_force>();
		for(int i=1;i<=A.n;++i)vs[i]=vec3{0,0,0};
		for(int i=1;i<=A.n;++i)if((A.info[i]&15)==2){
			int id=A.info[i]>>8;
			vs[i]=ri[id].v+FR[id]*(dt/ri[id].mass)+cross(
				TR[id]+ri[id].Iinv*(TR[id]+cross(ri[id].I*ri[id].w,ri[id].w))*dt,A.x[i]-ri[id].c
			);
		}
		for_neighbor<update_s>();
// if(stepp==258){
// 	int i;
// 	for(i=1;i<=A.n;++i)if(in_contact[i]){
// 		lp[i]=1;
// 		break;
// 	}
// 	++i;
// 	for(;i<=A.n;++i)lp[i]=0;
// }
		for_neighbor<update_rigid_pressure_force>();
		for(int i=1;i<=rcnt;++i)FR[i]=TR[i]=vec3{0,0,0};
		for(int i=1;i<=A.n;++i)if(in_contact[i]){
			update_rigid(i,tmp_ap[i]*A.mass[i]);
		}
		for(int i=1;i<=A.n;++i)if(in_contact[i]){
			break;
		}
		for(int i=1;i<=A.n;++i)if((A.info[i]&15)==2){
			int id=A.info[i]>>8;
			vr[i]=FR[id]*(dt/ri[id].mass)+cross(
				ri[id].Iinv*TR[id]*dt,A.x[i]-ri[id].c
			);
		}
		for_neighbor<update4>();
		memcpy(lp+1,p+1,A.n*4);
		sum_err0/=tot0;
		sum_err1/=tot1;
// cerr<<sum_err0<<' '<<sum_err1<<endl;
		if(tot0==0 || sum_err0<1e-2)flag=1;
		if(tot1==0 || sum_err1<0.1)break;
		if(ite>0 && last<=sum_err1)break;
		last=sum_err1;
// cerr<<"gua\n";
// cerr<<flag<<endl;
		sum_err0=0;
		tot0=0;
		sum_err1=0;
		tot1=0;
	}
	memcpy(FR+1,FR_nonp+1,rcnt*sizeof(vec3));
	memcpy(TR+1,TR_nonp+1,rcnt*sizeof(vec3));
	for_neighbor<update_fluid_pressure_force>();
	for(int i=1;i<=A.n;++i)if((A.info[i]&15)==2){
		int id=A.info[i]>>8;
		vs[i]=ri[id].v+FR[id]*(dt/ri[id].mass)+cross(
			TR[id]+ri[id].Iinv*(TR[id]+cross(ri[id].I*ri[id].w,ri[id].w))*dt,A.x[i]-ri[id].c
		);
	}
	for_neighbor<update_rigid_pressure_force>();
	for(int i=1;i<=A.n;++i)if(in_contact[i]){
		if(isnan(tmp_ap[i].x))
			++i,--i;
		update_rigid(i,tmp_ap[i]*A.mass[i]);
		if(fabs(FR[A.info[i]>>8].z)>100)
			++i,--i;
	}
	for(int i=1;i<=rcnt;++i){
		assert(!isnan(FR[i].x));
		ri[i].v=ri[i].v+FR[i]*(dt/ri[i].mass);
		ri[i].w=ri[i].w+ri[i].Iinv*(TR[i]+cross(ri[i].I*ri[i].w,ri[i].w))*dt;
		float nz=(ri[i].c+ri[i].v*dt).z;
		if(nz<h+0.25 && i==7)ri[i].v=ri[i].v+vec3{0,0,(h+0.25-nz)/dt};
		if(len(ri[i].w)>10)ri[i].w=ri[i].w*(10/len(ri[i].w));
		if(len(ri[i].v)>3)ri[i].v=ri[i].v*(3/len(ri[i].v));
		ri[i].trans=(I3+a2m(ri[i].w)*dt)*ri[i].trans;
	}
	for(int i=1;i<=A.n;++i)if((A.info[i]&15)!=1){
		float old_d=len(A.x[i]-ri[A.info[i]>>8].c);
		if((A.info[i]&15)==0){
			A.v[i]=A.v[i]+(fnonp[i]+tmp_ap[i])*(dt/A.mass[i]);
if(len(A.v[i])>20)A.v[i]=A.v[i]*(20/len(A.v[i]));
		}else{
			int id=A.info[i]>>8;
			A.v[i]=ri[id].v+cross(ri[id].w,A.x[i]-ri[id].c);
		}
if(isnan(A.v[i].x) || isnan(A.v[i].y) || isnan(A.v[i].z)){
	cerr<<A.v[i].x<<' '<<A.v[i].y<<' '<<A.v[i].z<<endl;
	exit(0);
}
if(isnan(A.x[i].x) || isnan(A.x[i].y) || isnan(A.x[i].z)){
	cerr<<A.x[i].x<<' '<<A.x[i].y<<' '<<A.x[i].z<<endl;
	exit(0);
}
//		cerr<<A.v[i].z<<endl;
//		cerr<<A.v[i].z<<endl;
		A.x[i]=A.x[i]+A.v[i]*dt;
		float new_d=len(A.x[i]-(ri[A.info[i]>>8].c+ri[A.info[i]>>8].v*dt));
		if(fabs(old_d-new_d)>1e-3){
			++i,--i;
		}
		if((A.info[i]&15)==0){
			for(int j=0;j<3;++j){				
				if(A.x[i][j]<lb[j]+h*2){
					A.x[i][j]=lb[j]+h*2;
					A.v[i][j]*=-0.5;
				}
				if(A.x[i][j]>rb[j]-h*2){
					A.x[i][j]=rb[j]-h*2;
					A.v[i][j]*=-0.5;
				}
			}
		}
// if((A.x[i]+A.v[i]*dt).x<lb.x || (A.x[i]+A.v[i]*dt).x>rb.x){
// 	cerr<<A.v[i].x<<' '<<A.x[i].x<<' '<<' '<<A.x[i].y<<' '<<A.x[i].z<<' '<<A.info[i]<<' '<<i<<endl;
// 	exit(0);
// }
// if((A.x[i]+A.v[i]*dt).y<lb.y || (A.x[i]+A.v[i]*dt).y>rb.y){
// 	cerr<<A.v[i].y<<' '<<A.x[i].x<<' '<<' '<<A.x[i].y<<' '<<A.x[i].z<<' '<<A.info[i]<<' '<<i<<endl;
// 	exit(0);
// }
if(A.x[i].z<lb.z || A.x[i].z>rb.z){
	cerr<<A.v[i].z<<' '<<A.x[i].x<<' '<<' '<<A.x[i].y<<' '<<A.x[i].z<<' '<<A.info[i]<<' '<<i<<endl;
	exit(0);
}
	}
	for(int i=1;i<=rcnt;++i)ri[i].c=ri[i].c+ri[i].v*dt;
}

char obuf[1<<24],*oh;
inline void write(float x){
	if(x<0){
		putchar('-');
		x=-x;
	}
	int t=floor(x);
	assert(t<10);
	*oh++=t+'0';
	*oh++='.';
	x-=t;
	for(int i=0;i<5;++i){
		x*=10;
		t=floor(x);
		*oh++=t+'0';
		x-=t;
	}
	*oh++=' ';
}
inline void out(vec3 P,vec3 C){
	write(P.x);
	write(P.z);
	write(P.y);
	*oh++=' ';
	write(C.x);
	write(C.y);
	write(C.z);
	*oh++='\n';
}

int main(){
	static int time_p[N],time_r[N2];
	freopen("scene.txt","r",stdin);
	int PC,RC,FC;
	scanf("%d",&PC);
	for(int i=1;i<=PC;++i)scanf("%f%f%f%f%f%f%f%f%d%d"
		,&A.x[i].x,&A.x[i].y,&A.x[i].z
		,&A.v[i].x,&A.v[i].y,&A.v[i].z
		,&A.mass[i],&A.rho[i],&A.info[i],&time_p[i]
	);
	scanf("%d",&RC);
	for(int i=1;i<=RC;++i)ri[i].input(),scanf("%d",&time_r[i]);

	float tt;
	float rt=0;
	scanf("%d%f",&FC,&tt);
	for(int fr=0;fr<FC;){
		rt-=dt;
		if(rt<=0){
			for(;A.n<PC && time_p[A.n+1]<=fr;++A.n);
			for(;rcnt<RC && time_r[rcnt+1]<=fr;++rcnt);
			rt=tt;
			char c[99];
			sprintf(c,"raw/%d.txt",fr);
			FILE*f=fopen(c,"w");
			oh=obuf;
			for(int i=1;i<=A.n;++i)if(A.info[i]!=1 || A.v[i].x==114514){
				out(A.x[i],colors[A.info[i]>>4&15]);
			}
			// for(int i=0;i<8;++i){
			// 	out(vec3{float(i>>2&1)*3,float(i>>1&1)*2,float(i&1)*2},WHITE);
			// }
			fwrite(obuf,1,oh-obuf,f);
			fclose(f);
			sprintf(c,"raw/%d_rigid.txt",fr);
			f=fopen(c,"w");
			fprintf(f,"%d\n",rcnt);
			for(int i=1;i<=rcnt;++i){
				for(int j=0;j<3;++j)fprintf(f,"%.5f%c",ri[i].c[j],j==2?'\n':' ');
				for(int j=0;j<3;++j)for(int k=0;k<3;++k)fprintf(f,"%.5f%c",ri[i].trans.a[j][k],k==2?'\n':' ');
			}
			fclose(f);
			++fr;
			fprintf(stderr,"frame:%d time:%d\n",fr,clock());
// break;
		}
		step();
	}
}