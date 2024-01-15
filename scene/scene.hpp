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

int time_p[N],time_r[N2];
void deal_rigid(int id,int be,int en){
	char c[99];
	sprintf(c,"tmp2/tmpp%d.txt",id);
	freopen(c,"w",stdout);
	for(int i=be;i<=en;++i)printf("%.6f %.6f %.6f\n",A.x[i].x,A.x[i].y,A.x[i].z);
}
void add(int type,int nx,int ny,int nz,float xs,float xl,float ys,float yl,float zs,float zl){
	if(type==2){
		++rcnt;
		ri[rcnt].mass=xl*yl*zl*(rho0/10);
		ri[rcnt].c=vec3{xs+xl/2,ys+yl/2,zs+zl/2};
		ri[rcnt].v=vec3{0,0,0};
		ri[rcnt].w=vec3{0,0,0};
		ri[rcnt].I=mat3{
			ri[rcnt].mass*(yl*yl+zl*zl)/12,0,0,
			0,ri[rcnt].mass*(xl*xl+zl*zl)/12,0,
			0,0,ri[rcnt].mass*(xl*xl+yl*yl)/12
		};
		ri[rcnt].Iinv=getinv(ri[rcnt].I);
		ri[rcnt].trans=I3;
	}
	float r=type!=1?h/2:h*2.5;
	// r=h/2;
	int oldn=A.n;
	for(int i=1;i<=nx;++i)
		for(int j=1;j<=ny;++j)
			for(int k=1;k<=nz;++k){
				++A.n;
				A.x[A.n]=vec3{xs+xl*i/(nx+1),ys+yl*j/(ny+1),zs+zl*k/(nz+1)};
				A.v[A.n]=vec3{0,0,0};
				A.mass[A.n]=pi*4/3*r*r*r*rho0;
				A.rho[A.n]=rho0;
				A.info[A.n]=type;
				if(type==2)A.info[A.n]+=(rcnt<<8)+(1<<4);
	}
	if(type==2){
		deal_rigid(rcnt,oldn+1,A.n);
	}
}
void add11(vec3 o,vec3 p1,vec3 p2,int color=1,vector<int>ban={}){
	p2=p2*(1.0/len(p2));
	vec3 p3=cross(p1,p2);
	float rrr=0.08;
	// for(int i=-2;i<=2;++i)for(int j=-2;j<=2;++j)if(abs(i)==2 || abs(j)==2){
	// 	++A.n;
	// 	A.x[A.n]=o+p1*i*d+p2*j*d;
	// 	A.v[A.n]=vec3{0,0,0};
	// 	A.rho[A.n]=rho0;
	// 	A.info[A.n]=2+(rcnt<<8)+(1<<4);
	// }
	for(int j:{-1,0,1})for(int k=0;k<3;++k)for(int i=0;i<48;++i){
		if(find(ban.begin(),ban.end(),i)!=ban.end())continue;
		++A.n;
		float alpha=2*pi*i/48;
		float r0=rrr+k*0.01;
		A.x[A.n]=o+p1*r0*cos(alpha)+p2*r0*sin(alpha)+p3*j*0.01;
		if(A.x[A.n].z<-h){
			--A.n;
			continue;
		}
		A.v[A.n]=vec3{0,0,0};
		A.rho[A.n]=rho0;
		A.info[A.n]=2+(rcnt<<8)+(color<<4);
	}
}
void add2(vec3 o,float R,vector<vec3>k){
	++rcnt;
	ri[rcnt].mass=pi*4/3*R*R*R*(rho0/10);
	ri[rcnt].c=o;
	ri[rcnt].v=vec3{0,0,0};
	ri[rcnt].w=vec3{0,0,0};
	ri[rcnt].I=I3*(ri[rcnt].mass*R*R*0.4);
	ri[rcnt].Iinv=getinv(ri[rcnt].I);
	ri[rcnt].trans=I3;
	float r=0.03;
	int t=ceil(R/r);
	int oldn=A.n;
	for(int i=-t;i<=t;++i)for(int j=-t;j<=t;++j)for(int k=-t;k<=t;++k){
		vec3 p=vec3{i,j,k}*r;
		if(len(p)>R)continue;
		if(len(p)<R-h*2)continue;
		++A.n;
		A.x[A.n]=o+p;
		A.v[A.n]=vec3{0,0,0};
		A.rho[A.n]=rho0;
		A.info[A.n]=2+(rcnt<<8)+(1<<4);
	}
	for(auto u:k)add11(o+u*(R+0.05),vec3{0,0,1},u);
	deal_rigid(rcnt,oldn+1,A.n);
}
void add3(vec3 o,vec3 p1,vec3 p2,int color=1,vector<int>ban={}){
	++rcnt;
	int oldn=A.n;
	add11(o,p1,p2,color,ban);
	ri[rcnt].c=o;
	ri[rcnt].v=vec3{0,0,0};
	ri[rcnt].w=vec3{0,0,0};
	for(int i=oldn+1;i<=A.n;++i){
		vec3 p=A.x[i]-o;
		mat3 dI=mat3{
			sqr(p.y)+sqr(p.z),-p.x*p.y,-p.x*p.z,
			-p.x*p.y,sqr(p.x)+sqr(p.z),-p.y*p.z,
			-p.x*p.z,-p.y*p.z,sqr(p.x)+sqr(p.y)
		};
		ri[rcnt].I=ri[rcnt].I+dI;
	}
	ri[rcnt].mass=rho0*0.1*pi*h*(sqr(0.11)-sqr(0.06));
	ri[rcnt].I=ri[rcnt].I*ri[rcnt].mass;
	ri[rcnt].Iinv=getinv(ri[rcnt].I);
	ri[rcnt].trans=I3;
	deal_rigid(rcnt,oldn+1,A.n);
}