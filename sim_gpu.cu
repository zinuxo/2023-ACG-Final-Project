#include<cmath>
#include<iostream>
#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<algorithm>
#include<ctime>
#include<random>
#include<cassert>
#include <type_traits>  
#include<functional>
using namespace std;

constexpr int N=2e5+5,N2=105;

constexpr float pi=3.14159265358979;
constexpr float eps=1e-4;
constexpr float g=10;
constexpr float mu=0.003;
constexpr float rho0=1e3;
constexpr float k=40;
constexpr float dt=1e-3;
// constexpr float dt=4e-4;
constexpr float h=0.04;
constexpr float omega=0.05;
constexpr float gamma=0.7;


__device__ __host__ inline float sqr(float x){return x*x;}
__device__ __host__ inline float cubic(float x){return x*x*x;}
struct vec3{
	float x,y,z;
	inline float&operator[](int i){return i==0?x:i==1?y:z;}
	inline float operator[](int i)const{return i==0?x:i==1?y:z;}
};
__device__ __host__ inline float dot(vec3 a,vec3 b){return a.x*b.x+a.y*b.y+a.z*b.z;}
__device__ __host__ inline vec3 operator-(vec3 a,vec3 b){return vec3{a.x-b.x,a.y-b.y,a.z-b.z};}
__device__ __host__ inline bool operator==(vec3 a,vec3 b){return a.x==b.x && a.y==b.y && a.z==b.z;}
__device__ __host__ inline vec3 operator+(vec3 a,vec3 b){return vec3{a.x+b.x,a.y+b.y,a.z+b.z};}
__device__ __host__ inline vec3 operator*(vec3 a,float b){return vec3{a.x*b,a.y*b,a.z*b};}
__device__ __host__ inline float len2(vec3 a){return dot(a,a);}
__device__ __host__ inline float len(vec3 a){return sqrt(len2(a));}
__device__ __host__ inline vec3 norm(vec3 a){return a*(1.0/len(a));}
__device__ __host__ inline vec3 cross(vec3 a,vec3 b){return vec3{a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
__device__ __host__ inline float getW(vec3 r){
	float q=len(r)/h;
	if(q<=1){
		return 3/(2*pi*h*h*h)*(2.0/3-q*q+0.5*q*q*q);
	}else if(q<=2){
		return 3/(2*pi*h*h*h)*1.0/6*(2-q)*(2-q)*(2-q);
	}else return 0;
}
__device__ __host__ inline vec3 getWG(vec3 r){
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
constexpr vec3 lb=vec3{-h,-h,-h};
constexpr vec3 rb=vec3{3+h,2+h,2+h};
constexpr int M1=(rb.x-lb.x)/h/2+1,M2=(rb.y-lb.y)/h/2+1,M3=(rb.z-lb.z)/h/2+1,M=M1*M2*M3;

struct mat3{
	float a[3][3];
};
__device__ __host__ inline mat3 operator*(mat3 a,mat3 b){
	mat3 c;
	for(int i=0;i<3;++i)for(int j=0;j<3;++j){
		c.a[i][j]=0;
		for(int k=0;k<3;++k){
			c.a[i][j]+=a.a[i][k]*b.a[k][j];
		}
	}
	return c;
}
__device__ __host__ inline vec3 operator*(mat3 a,vec3 b){
	vec3 c;
	c.x = a.a[0][0] * b.x + a.a[0][1] * b.y + a.a[0][2] * b.z;
	c.y = a.a[1][0] * b.x + a.a[1][1] * b.y + a.a[1][2] * b.z;
	c.z = a.a[2][0] * b.x + a.a[2][1] * b.y + a.a[2][2] * b.z;
	return c;
}
__device__ __host__ inline mat3 operator+(mat3 a,mat3 b){
	mat3 c;
	for(int i=0;i<3;++i)for(int j=0;j<3;++j){
		c.a[i][j]=a.a[i][j]+b.a[i][j];
	}
	return c;
}
__device__ __host__ inline mat3 operator*(mat3 a,float b){
	mat3 c;
	for(int i=0;i<3;++i)for(int j=0;j<3;++j){
		c.a[i][j]=a.a[i][j]*b;
	}
	return c;
}
__device__ __host__ inline mat3 operator-(mat3 a,mat3 b){
	mat3 c;
	for(int i=0;i<3;++i)for(int j=0;j<3;++j){
		c.a[i][j]=a.a[i][j]-b.a[i][j];
	}
	return c;
}
__device__ __host__ inline mat3 getinv(mat3 a){
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
__device__ __host__ inline mat3 a2m(vec3 a){
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
constexpr mat3 I3=mat3{1,0,0,0,1,0,0,0,1};


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
	void move1(){
		memmove(x,x+1,n*sizeof(vec3));
		memmove(v,v+1,n*sizeof(vec3));
		memmove(mass,mass+1,n*sizeof(float));
		memmove(rho,rho+1,n*sizeof(float));
		memmove(info,info+1,n*sizeof(int));
	}
	void move2(){
		memmove(x+1,x,n*sizeof(vec3));
		memmove(v+1,v,n*sizeof(vec3));
		memmove(mass+1,mass,n*sizeof(float));
		memmove(rho+1,rho,n*sizeof(float));
		memmove(info+1,info,n*sizeof(int));
	}
}A,B;
__host__ __device__ __device__ __host__ inline int get_id(vec3 p){
	int x=floor((p.x-lb.x)/h/2);
	int y=floor((p.y-lb.y)/h/2);
	int z=floor((p.z-lb.z)/h/2);
	// if(0<=x && x<M1);else printf("%.6f %.6f %.6f\n",p.x,p.y,p.z);
	assert(0<=x && x<M1);
	assert(0<=y && y<M2);
	assert(0<=z && z<M3);
	return x*M2*M3+y*M3+z;
}

namespace CUDA{
__device__ Particles A,B;
__device__ Rigid ri[N2];
__device__ int rcnt;
typedef unsigned long long u64;
constexpr int LB=6,HB=10,U=1<<LB+HB;
__device__ u64 buffer_a[N],buffer_b[N];
__device__ int bucket[1<<HB],buffer_tmp[N],buc_be[M],buc_num[M];
__global__ void clear_buc(){
	int i1=threadIdx.y*blockDim.x+threadIdx.x;
	int i2=blockIdx.y*gridDim.x+blockIdx.x;
	int i=i2*blockDim.x*blockDim.y+i1;
	if(i<1<<HB)bucket[i]=0;
}
__global__ void mysort1(){
	int i1=threadIdx.y*blockDim.x+threadIdx.x;
	int i2=blockIdx.y*gridDim.x+blockIdx.x;
	int i=i2*blockDim.x*blockDim.y+i1;
	if(i<M)buc_num[i]=0;
	if(i>=A.n)return;
	u64 _val=get_id(A.x[i])|u64(i)<<16;
	buffer_a[i]=_val;
	int val=_val&(U-1);
	buffer_tmp[i]=atomicAdd(bucket+(val>>LB),1);
}
__global__ void mysort2(){
	int i1=threadIdx.y*blockDim.x+threadIdx.x;
	int i2=blockIdx.y*gridDim.x+blockIdx.x;
	int i=i2*blockDim.x*blockDim.y+i1;
	++i;
	if(i>1<<HB)return;
	__shared__ int c[1<<HB];
	c[i-1]=bucket[i-1];
	__syncthreads();
	for(int j=1;j<1<<HB;j<<=1){
		if((i&-i)==j)c[i-1+j]+=c[i-1];
		__syncthreads();
	}
	for(int j=1<<HB-1;j>=1;j>>=1){
		if((i&-i)==j && i>j)c[i-1]+=c[i-1-j];
		__syncthreads();
	}
	bucket[i-1]=c[i-1];
}
__global__ void mysort3(){
	int i1=threadIdx.y*blockDim.x+threadIdx.x;
	int i2=blockIdx.y*gridDim.x+blockIdx.x;
	int i=i2*blockDim.x*blockDim.y+i1;
	if(i>=A.n)return;
	int val=buffer_a[i]&(U-1);
	int hb=val>>LB;
	int offset=hb==0?0:bucket[hb-1];
	buffer_b[offset+buffer_tmp[i]]=buffer_a[i];
}
__global__ void mysort4(){
	int i1=threadIdx.y*blockDim.x+threadIdx.x;
	int i2=blockIdx.y*gridDim.x+blockIdx.x;
	if(i2>=1<<HB)return;
	int be=i2==0?0:bucket[i2-1],en=bucket[i2]-1;
	if(be+i1>en)return;
	static constexpr int N=1000;
	assert(en-be+1<=N);
	assert(en-be+1<=blockDim.x*blockDim.y);
	int num=en-be+1;
	int now=0;
	__shared__ u64 ta[2][N];
	__shared__ short cc[N][4];
	ta[now][i1]=buffer_b[be+i1];
	for(int bit=0;bit<LB;bit+=2){
		for(int i=0;i<4;++i)cc[i1][i]=0;
		int index=ta[now][i1]>>bit&3;
		cc[i1][index]=1;
		__syncthreads();
		int i=i1+1,j;
		for(j=1;j<num;j<<=1){
			if((i&-i)==j && i+j<=num)for(int k=0;k<4;++k)cc[i-1+j][k]+=cc[i-1][k];
			__syncthreads();
		}
		for(j>>=1;j>=1;j>>=1){
			if((i&-i)==j && i>j)for(int k=0;k<4;++k)cc[i-1][k]+=cc[i-1-j][k];
			__syncthreads();
		}
		int position=cc[i1][index];
		for(int i=0;i<index;++i)position+=cc[num-1][i];
		ta[now^1][position-1]=ta[now][i1];
		__syncthreads();
		now^=1;
	}
	buffer_a[be+i1]=ta[now][i1];
	if(i1==0 || (ta[now][i1]&(U-1))!=(ta[now][i1-1]&(U-1))){
		buc_be[ta[now][i1]&(U-1)]=be+i1;
		int cnt=1;
		for(;i1+cnt<num && (ta[now][i1+cnt]&(U-1))==(ta[now][i1]&(U-1));++cnt);
		buc_num[ta[now][i1]&(U-1)]=cnt;
	}
}
__global__ void a_to_b(){
	int i1=threadIdx.y*blockDim.x+threadIdx.x;
	int i2=blockIdx.y*gridDim.x+blockIdx.x;
	int i=i2*blockDim.x*blockDim.y+i1;
	if(i>=A.n)return;
	int id=buffer_a[i]>>16;
	B.x[i]=A.x[id];
	B.v[i]=A.v[id];
	B.mass[i]=A.mass[id];
	B.rho[i]=A.rho[id];
	B.info[i]=A.info[id];
}
__global__ void b_to_a(){
	int i1=threadIdx.y*blockDim.x+threadIdx.x;
	int i2=blockIdx.y*gridDim.x+blockIdx.x;
	int i=i2*blockDim.x*blockDim.y+i1;
	if(i>=A.n)return;
	A.x[i]=B.x[i];
	A.v[i]=B.v[i];
	A.mass[i]=B.mass[i];
	A.rho[i]=B.rho[i];
	A.info[i]=B.info[i];
	assert(get_id(B.x[i-1])<=get_id(B.x[i]));
}
inline void radix_sort(){
	static const dim3 grid(32,32),block(32,32),block_buc(32,32);
	clear_buc<<<dim3(1),block_buc>>>();
	mysort1<<<grid,block>>>();
	mysort2<<<dim3(1),block_buc>>>();
	mysort3<<<grid,block>>>();
	mysort4<<<grid,block>>>();
	a_to_b<<<grid,block>>>();
	b_to_a<<<grid,block>>>();
}


__device__ float aii[N],s[N],p[N],lp[N];
__device__ vec3 v2[N],fnonp[N],tmp_ap[N];
__device__ vec3 FR_nonp[N2],TR_nonp[N2],FR[N2],TR[N2],vs[N],vr[N],sumf[N];
__device__ bool has_contact[N],has_fluid[N];
__device__ int tot0,tot1;
__device__ float sum_err0,sum_err1;
__device__ inline void atomicAddvec3(vec3&a,vec3 b){
	atomicAdd(&a.x,b.x);
	atomicAdd(&a.y,b.y);
	atomicAdd(&a.z,b.z);
}
__device__ inline void update_rigid(int i,vec3 f){
	int id=A.info[i]>>8;--id;
	atomicAddvec3(FR[id],f);
	atomicAddvec3(TR[id],cross(A.x[i]-ri[id].c,f));
}
constexpr int MX=M1,MY=M2,MZ=M3;
__device__ inline mat3 getK(int i,int j){
	int R=A.info[i]>>8;--R;
	return I3*(1.0/ri[R].mass)-a2m(A.x[i]-ri[R].c)*ri[R].Iinv*a2m(A.x[j]-ri[R].c);
}
__device__ inline bool check(int type,int x,int y,int z){
	if(x<0 || x>=MX || y<0 || y>=MY || z<0 || z>=MZ)return false;
	// int id=z*MY*MX+y*MX+x;
	int id=x*M2*M3+y*M3+z;
	if(type==0)return 1;
	if(type==1)return has_fluid[id] || has_contact[id];
	if(type==2)return has_fluid[id];
	if(type==3)return has_contact[id];
	return false;
}
__device__ inline bool check2(int type,int info){
	if(type==0)return (info&15)!=1;
	if(type==1)return (info&15)==0 || (info&128);
	if(type==2)return (info&15)==0;
	if(type==3)return (info&128);
	return false;
}
template<typename T>constexpr __device__ float*get_fp(int i){
	return T::fp[i];
}
template<typename T>constexpr __device__ vec3*get_vp(int i){
	return T::vp[i];
}
template<typename T>constexpr __device__ int*get_ip(int i){
	return T::ip[i];
}
template<
	typename T,
	void (*init)(T&,int,float**,vec3**,int**),
	void (*update)(T&,int,int),
	void (*finish)(T&,int,int)
>__global__ void for_neighbor(){
	static constexpr int N=500;
	__shared__ float sf[3][N];
	__shared__ vec3 sv[3][N];
	__shared__ int si[1][N];
	int i1=threadIdx.y*blockDim.x+threadIdx.x;
	int y=blockIdx.y,z=blockIdx.z;
	for(int x=0;x<MX;++x){
		if(!check(T::check_type,x,y,z))continue;
		int myid=-1,myid2,mynum=0,num=0;
		for(int dx=-1;dx<=1;++dx)for(int dy=-1;dy<=1;++dy)for(int dz=-1;dz<=1;++dz){
			int xx=x+dx,yy=y+dy,zz=z+dz;
			if(xx<0 || xx>=MX || yy<0 || yy>=MY || zz<0 || zz>=MZ)continue;
			if(!check(T::check_type,xx,yy,zz))continue;
			// int id=zz*MY*MX+yy*MX+xx;
			int id=xx*M2*M3+yy*M3+zz;
			int _be=buc_be[id],_num=buc_num[id];
			if(dx==0 && dy==0 && dz==0)myid=_be,mynum=_num,myid2=num;
// if(num+_num>500)printf("%d %d\n",num,_num);
			assert(num+_num<=N);
			assert(_num<=blockDim.x*blockDim.y);
			if(i1<_num){
				sf[0][num+i1]=get_fp<T>(0)[_be+i1];
				if(0<T::fpl)sf[0][num+i1]=get_fp<T>(0)[_be+i1];
				if(1<T::fpl)sf[1][num+i1]=get_fp<T>(1)[_be+i1];
				if(2<T::fpl)sf[2][num+i1]=get_fp<T>(2)[_be+i1];
				if(0<T::vpl)sv[0][num+i1]=get_vp<T>(0)[_be+i1];
				if(1<T::vpl)sv[1][num+i1]=get_vp<T>(1)[_be+i1];
				if(2<T::vpl)sv[2][num+i1]=get_vp<T>(2)[_be+i1];
				if(0<T::ipl)si[0][num+i1]=get_ip<T>(0)[_be+i1];
			}
			num+=_num;
		}
		// int load_id=-1;
		// for(int dx=-1;dx<=1;++dx)for(int dy=-1;dy<=1;++dy)for(int dz=-1;dz<=1;++dz){
		// 	int xx=x+dx,yy=y+dy,zz=z+dz;
		// 	if(xx<0 || xx>=MX || yy<0 || yy>=MY || zz<0 || zz>=MZ)continue;
		// 	if(!check(T::check_type,xx,yy,zz))continue;
		// 	int id=xx*M2*M3+yy*M3+zz;
		// 	int _be=buc_be[id],_num=buc_num[id];
		// 	if(dx==0 && dy==0 && dz==0)myid=_be,mynum=_num,myid2=num;
		// 	assert(num+_num<=N);
		// 	assert(_num<=blockDim.x*blockDim.y);
		// 	if(num<=i1 && i1<num+_num)load_id=_be+i1-num;
		// 	num+=_num;
		// }
		// __syncthreads();
		// assert(num<=blockDim.x*blockDim.y);
		// assert(myid!=-1);
		// if(load_id!=-1){
		// 	assert(i1<N);
		// 	assert(load_id<A.n);
		// 	sf[0][i1]=get_fp<T>(0)[load_id];
		// 	if(0<T::fpl)sf[0][i1]=get_fp<T>(0)[load_id];
		// 	if(1<T::fpl)sf[1][i1]=get_fp<T>(1)[load_id];
		// 	if(2<T::fpl)sf[2][i1]=get_fp<T>(2)[load_id];
		// 	if(0<T::vpl)sv[0][i1]=get_vp<T>(0)[load_id];
		// 	if(1<T::vpl)sv[1][i1]=get_vp<T>(1)[load_id];
		// 	if(2<T::vpl)sv[2][i1]=get_vp<T>(2)[load_id];
		// 	if(0<T::ipl)si[0][i1]=get_ip<T>(0)[load_id];
		// }
		__syncthreads();
		if(i1<mynum && check2(T::check_type,si[0][myid2+i1])){
			T local_updater;
			float*fp[3]={sf[0],sf[1],sf[2]};
			vec3*vp[3]={sv[0],sv[1],sv[2]};
			int*ip[1]={si[0]};
			init(local_updater,i1,fp,vp,ip);
			for(int j=0;j<num;++j)if(len2(sv[0][myid2+i1]-sv[0][j])<=sqr(2*h)){
				update(local_updater,myid2+i1,j);
			}
			finish(local_updater,myid2+i1,myid+i1);
		}
	}
}
struct update_rho{
	static constexpr int check_type=0;
	static constexpr int fpl=1,vpl=1,ipl=1;
	static constexpr float*fp[3]={A.mass,NULL,NULL};
	static constexpr vec3*vp[3]={A.x,NULL,NULL};
	static constexpr int*ip[1]={A.info};
	float*mass;
	vec3*x;
	int*info;
	float new_rho;
	bool in_contact;
};
__device__ inline void update_rho_init(update_rho&A,int i,float**fp,vec3**vp,int**ip){
	A.new_rho=0;
	A.in_contact=false;
	A.mass=fp[0];
	A.x=vp[0];
	A.info=ip[0];
}
__device__ inline void update_rho_update(update_rho&A,int i,int j){
	A.new_rho+=A.mass[j]*getW(A.x[i]-A.x[j]);
	if((A.info[i]&15)!=0 && (A.info[j]&15)!=0 && (A.info[i]>>8)!=(A.info[j]>>8))
		A.in_contact=true;
}
__device__ inline void update_rho_finish(update_rho&U,int i,int ii){
	U.new_rho=max(U.new_rho,rho0);
	A.rho[ii]=U.new_rho;
	int id=get_id(U.x[i]);
	if(U.in_contact){
		A.info[ii]|=128;
		has_contact[id]=true;
	}
	if((U.info[i]&15)==0){
		has_fluid[id]=true;
	}
}
struct update2{
	static constexpr int check_type=1;
	static constexpr int fpl=2,vpl=2,ipl=1;
	static constexpr float*fp[3]={A.mass,A.rho,NULL};
	static constexpr vec3*vp[3]={A.x,A.v,NULL};
	static constexpr int*ip[1]={A.info};
	float*mass,*rho;
	vec3*x,*v;
	int*info;
	vec3 sum1,nonp;
	float sum2;
	vec3 sum_f;
};
__device__ inline void update2_init(update2&U,int i,float**fp,vec3**vp,int**ip){
	U.sum1=vec3{0,0,0};
	U.sum2=0;
	U.sum_f=vec3{0,0,0};
	U.mass=fp[0];
	U.rho=fp[1];
	U.x=vp[0];
	U.v=vp[1];
	U.info=ip[0];

	U.nonp=vec3{0,0,-g*U.mass[i]};
}
__device__ inline void update2_update(update2&A,int i,int j){
	if(len(A.x[i]-A.x[j])<eps)return;
	if((A.info[i]&15)==0){
		A.sum1=A.sum1+getWG(A.x[i]-A.x[j])*A.mass[j];
		A.sum2=A.sum2+len2(getWG(A.x[i]-A.x[j]))*A.mass[j];
		vec3 f=(A.v[i]-A.v[j])*(mu*(A.mass[j]/A.rho[j])*2*
			len(getWG(A.x[i]-A.x[j]))/len(A.x[i]-A.x[j]));
		A.nonp=A.nonp-f;
		if((A.info[j]&15)==2){
			int id=A.info[j]>>8;--id;
			atomicAddvec3(FR[id],f);
			atomicAddvec3(TR[id],cross(A.x[j]-ri[id].c,f));
		}
	}else if((A.info[j]&15)!=0 && (A.info[i]>>8)!=(A.info[j]>>8)){
		A.sum_f=A.sum_f-getWG(A.x[i]-A.x[j])*(A.mass[j]*(1/sqr(A.rho[i]))*A.mass[i]*dt);
	}
}
__device__ inline void update2_finish(update2&U,int i,int ii){
	if((U.info[i]&15)==0){
		aii[ii]=-(len2(U.sum1)+U.sum2*U.mass[i])*sqr(dt)/sqr(U.rho[i]);
		fnonp[ii]=U.nonp;
		v2[ii]=U.v[i]+U.nonp*(dt/U.mass[i]);
	}else{
		sumf[ii]=U.sum_f;
	}
}
struct update3{
	static constexpr int check_type=1;
	static constexpr int fpl=1,vpl=3,ipl=1;
	static constexpr float*fp[3]={A.mass,NULL,NULL};
	static constexpr vec3*vp[3]={A.x,v2,sumf};
	static constexpr int*ip[1]={A.info};
	float*mass;
	vec3*x,*v2,*sumf;
	int*info;
	float sum;
};
__device__ inline void update3_init(update3&U,int i,float**fp,vec3**vp,int**ip){
	U.sum=0;
	U.mass=fp[0];
	U.x=vp[0];
	U.v2=vp[1];
	U.sumf=vp[2];
	U.info=ip[0];
}
__device__ inline void update3_update(update3&A,int i,int j){
	if((A.info[i]&15)==0){
		A.sum+=A.mass[j]*dot(A.v2[i]-A.v2[j],getWG(A.x[i]-A.x[j]));
	}else if(i!=j && (A.info[j]&15)!=0){
		A.sum-=dot(((A.info[i]>>8)==(A.info[j]>>8)?getK(j,i)*A.sumf[i]:vec3{0,0,0})
			-getK(i,i)*A.sumf[i],getWG(A.x[i]-A.x[j]))*A.mass[j];
	}
}
__device__ inline void update3_finish(update3&U,int i,int ii){
	if((U.info[i]&15)==0){
		s[ii]=rho0-A.rho[ii]-U.sum*dt;
	}else{
		aii[ii]=U.sum;
	}		
}
struct update_fluid_pressure_force{
	static constexpr int check_type=2;
	static constexpr int fpl=3,vpl=1,ipl=1;
	static constexpr float*fp[3]={A.mass,A.rho,lp};
	static constexpr vec3*vp[3]={A.x,NULL,NULL};
	static constexpr int*ip[1]={A.info};
	float*mass,*rho,*lp;
	vec3*x;
	int*info;
	vec3 sum;
};
__device__ inline void update_fluid_pressure_force_init(update_fluid_pressure_force&U,int i,float**fp,vec3**vp,int**ip){
	U.sum=vec3{0,0,0};
	U.mass=fp[0];
	U.rho=fp[1];
	U.lp=fp[2];
	U.x=vp[0];
	U.info=ip[0];
}
__device__ inline void update_fluid_pressure_force_update(update_fluid_pressure_force&A,int i,int j){
	vec3 f=getWG(A.x[i]-A.x[j])*(A.mass[j]*(A.lp[i]/sqr(A.rho[i])+
		((A.info[j]&15)==0?A.lp[j]/sqr(A.rho[j]):0)));
	A.sum=A.sum-f;
	if((A.info[j]&15)==2){
		int id=A.info[j]>>8;--id;
		atomicAddvec3(FR[id],f);
		atomicAddvec3(TR[id],cross(A.x[j]-ri[id].c,f));
	}
}
__device__ inline void update_fluid_pressure_force_finish(update_fluid_pressure_force&U,int i,int ii){
	tmp_ap[ii]=U.sum;
// if(len2(U.sum)>1e-6)assert(0);
}
struct update_s{
	static constexpr int check_type=3;
	static constexpr int fpl=1,vpl=2,ipl=1;
	static constexpr float*fp[3]={A.mass,NULL,NULL};
	static constexpr vec3*vp[3]={A.x,vs,NULL};
	static constexpr int*ip[1]={A.info};
	float*mass;
	vec3*x,*vs;
	int*info;
	float sum;
};
__device__ inline void update_s_init(update_s&U,int i,float**fp,vec3**vp,int**ip){
	U.sum=0;
	U.mass=fp[0];
	U.x=vp[0];
	U.vs=vp[1];
	U.info=ip[0];
}
__device__ inline void update_s_update(update_s&A,int i,int j){
	if(i!=j && (A.info[j]&15)!=1){
		A.sum=A.sum+dot(A.vs[j]-A.vs[i],getWG(A.x[i]-A.x[j]))*A.mass[j];
	}
}
__device__ inline void update_s_finish(update_s&U,int i,int ii){
	s[ii]=(rho0-A.rho[ii])/dt+U.sum;
}
struct update_rigid_pressure_force{
	static constexpr int check_type=3;
	static constexpr int fpl=3,vpl=1,ipl=1;
	static constexpr float*fp[3]={A.mass,A.rho,lp};
	static constexpr vec3*vp[3]={A.x,NULL,NULL};
	static constexpr int*ip[1]={A.info};
	float*mass,*rho,*lp;
	vec3*x;
	int*info;
	vec3 sum;
};
__device__ inline void update_rigid_pressure_force_init(update_rigid_pressure_force&U,int i,float**fp,vec3**vp,int**ip){
	U.sum=vec3{0,0,0};
	U.mass=fp[0];
	U.rho=fp[1];
	U.lp=fp[2];
	U.x=vp[0];
	U.info=ip[0];
}
__device__ inline void update_rigid_pressure_force_update(update_rigid_pressure_force&A,int i,int j){
	
	if((A.info[i]>>8)!=(A.info[j]>>8) && (A.info[j]&15)!=0){
		vec3 f=getWG(A.x[i]-A.x[j])*(A.mass[j]*(A.lp[i]/sqr(A.rho[i])+
			((A.info[j]&15)==2?A.lp[j]/sqr(A.rho[j]):0)));
		A.sum=A.sum-f;
	}
}
__device__ inline void update_rigid_pressure_force_finish(update_rigid_pressure_force&U,int i,int ii){
	
	tmp_ap[ii]=U.sum;
}
__device__ bool flag;
struct update4{
	static constexpr int check_type=1;
	static constexpr int fpl=1,vpl=3,ipl=1;
	static constexpr float*fp[3]={A.mass,NULL,NULL};
	static constexpr vec3*vp[3]={A.x,vr,tmp_ap};
	static constexpr int*ip[1]={A.info};
	float*mass;
	vec3*x,*vr,*tmp_ap;
	int*info;
	float sum=0;
};
__device__ inline void update4_init(update4&U,int i,float**fp,vec3**vp,int**ip){
	U.sum=0;
	U.mass=fp[0];
	U.x=vp[0];
	U.vr=vp[1];
	U.tmp_ap=vp[2];
	U.info=ip[0];
}
__device__ inline void update4_update(update4&A,int i,int j){
	if((A.info[i]&128)){
		if((A.info[j]&15)==0)return;
		A.sum=A.sum-dot(((A.info[i]>>8)==(A.info[j]>>8)?A.vr[j]:vec3{0,0,0})-A.vr[i],getWG(A.x[i]-A.x[j]))*A.mass[j];
	}else{
		A.sum=A.sum+A.mass[j]*dot(A.tmp_ap[i]-((A.info[j]&15)==0?A.tmp_ap[j]:vec3{0,0,0}),getWG(A.x[i]-A.x[j]));
	}
}
__device__ inline void update4_finish(update4&U,int i,int ii){
	
	float sum=U.sum;
	if((U.info[i]&15)==0){
		sum*=sqr(dt);
	}
	if((U.info[i]&128) && s[ii]==0)return;
// assert(flag);
	if((U.info[i]&15)==0 && flag){
		return;
	}
// printf("%d %d\n",U.info[i],ii);
// assert((U.info[i]&128)==0);
// assert(ii!=1047);
// assert(s[ii]==sum);
// assert(s[ii]==0);
// assert(sum==0);
// if(s[ii]!=sum)printf("%.6f %.6f\n",s[ii],sum);
	if(fabs(aii[ii])>((U.info[i]&128)?5e-3:1e-5)){
		p[ii]=max(lp[ii]+((U.info[i]&128)?0.1f:omega)/aii[ii]*(s[ii]-sum),0.0f);
	}else p[ii]=0;
// if(p[ii]!=0)assert(0);
	float err=fabs(s[ii]-sum)/rho0;
	float&sum_err=(U.info[i]&15)==0?sum_err0:sum_err1;
	int&tot=(U.info[i]&15)==0?tot0:tot1;
	sum_err+=err;
	++tot;
}

template<void (*upda)(int)>__global__ void _forP(){
	int i1=threadIdx.y*blockDim.x+threadIdx.x;
	int i2=blockIdx.y*gridDim.x+blockIdx.x;
	int i=i2*blockDim.x*blockDim.y+i1;
	if(i>=A.n)return;
	upda(i);
}
template<void (*upda)(int)>void forP(){
	static const dim3 grid=dim3(16,16),block=dim3(32,32);
	_forP<upda><<<grid,block>>>();
}
template<void (*upda)(int)>__global__ void _forR(){
	int i1=threadIdx.y*blockDim.x+threadIdx.x;
	int i2=blockIdx.y*gridDim.x+blockIdx.x;
	int i=i2*blockDim.x*blockDim.y+i1;
	if(i>=rcnt)return;
	upda(i);
}
template<void (*upda)(int)>void forR(){
	static const dim3 grid=dim3(1,1),block=dim3(8,8);
	_forR<upda><<<grid,block>>>();
}
__device__ void upda1(int i){
	A.info[i]&=~128;
	v2[i]=tmp_ap[i]=vr[i]=vec3{0,0,0};
	lp[i]=p[i]=aii[i]=0;
	has_contact[i]=has_fluid[i]=false;
}
__device__ void upda2(int i){	
	FR[i]=vec3{0,0,-g*ri[i].mass},TR[i]=vec3{0,0,0};
}
__device__ void upda3(int i){
	FR_nonp[i]=FR[i];
	TR_nonp[i]=TR[i];
}
__device__ void upda4(int i){
	FR[i]=FR_nonp[i];
	TR[i]=TR_nonp[i];
}
__device__ void upda5(int i){
	if((A.info[i]&15)==2){
		int id=A.info[i]>>8;--id;
		vs[i]=ri[id].v+FR[id]*(dt/ri[id].mass)+cross(
			TR[id]+ri[id].Iinv*(TR[id]+cross(ri[id].I*ri[id].w,ri[id].w))*dt,A.x[i]-ri[id].c
		);
	}else vs[i]=vec3{0,0,0};
}
__device__ void upda6(int i){
	FR[i]=TR[i]=vec3{0,0,0};
}
__device__ void upda7(int i){
	if((A.info[i]&128)){
		CUDA::update_rigid(i,tmp_ap[i]*A.mass[i]);
	}
}
__device__ void upda8(int i){
	if((A.info[i]&15)==2){
		int id=A.info[i]>>8;--id;
		vr[i]=FR[id]*(dt/ri[id].mass)+cross(
			ri[id].Iinv*TR[id]*dt,A.x[i]-ri[id].c
		);
	}
}
__device__ void upda9(int i){
	lp[i]=p[i];
}
__device__ void upda10(int i){
	ri[i].v=ri[i].v+FR[i]*(dt/ri[i].mass);
	ri[i].w=ri[i].w+ri[i].Iinv*(TR[i]+cross(ri[i].I*ri[i].w,ri[i].w))*dt;
	if(len(ri[i].w)>10)ri[i].w=ri[i].w*(10/len(ri[i].w));
	if(len(ri[i].v)>3)ri[i].v=ri[i].v*(3/len(ri[i].v));
}
__device__ void upda11(int i){
	if((A.info[i]&15)!=1){
// float old=A.v[i].z;
		if((A.info[i]&15)==0){
			A.v[i]=A.v[i]+(fnonp[i]+tmp_ap[i])*(dt/A.mass[i]);
		}else{
			int id=A.info[i]>>8;--id;
			A.v[i]=ri[id].v+cross(ri[id].w,A.x[i]-ri[id].c);
		}
// printf("%d %d %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n",i,A.info[i],old,A.v[i].z,fnonp[i].z,tmp_ap[i].z,dt/A.mass[i],(fnonp[i]+tmp_ap[i]).z*(dt/A.mass[i]),A.mass[i]);
		A.x[i]=A.x[i]+A.v[i]*dt;
		if((A.info[i]&15)==0){
			if(A.x[i].x<lb.x+h*2){
				A.x[i].x=lb.x+h*2;
				A.v[i].x*=-0.5;
			}
			if(A.x[i].x>rb.x-h*2){
				A.x[i].x=rb.x-h*2;
				A.v[i].x*=-0.5;
			}
			if(A.x[i].y<lb.y+h*2){
				A.x[i].y=lb.y+h*2;
				A.v[i].y*=-0.5;
			}
			if(A.x[i].y>rb.y-h*2){
				A.x[i].y=rb.y-h*2;
				A.v[i].y*=-0.5;
			}
			if(A.x[i].z<lb.z+h*2){
				A.x[i].z=lb.z+h*2;
				A.v[i].z*=-0.5;
			}
			if(A.x[i].z>rb.z-h*2){
				A.x[i].z=rb.z-h*2;
				A.v[i].z*=-0.5;
			}
		}
	}
}
__device__ void upda12(int i){
	ri[i].c=ri[i].c+ri[i].v*dt;
}

__device__ float last;
__global__ void init_ite(){
	
	sum_err0=0;
	tot0=0;
	sum_err1=0;
	tot1=0;
	flag=false;
}
__device__ bool to_stop;
bool to_stop_host;
__global__ void stop_ite(int ite){
	
	sum_err0/=tot0;
	sum_err1/=tot1;
	if(tot0==0 || sum_err0<1e-2)flag=1;
	to_stop=false;
	if(tot1==0 || sum_err1<0.1)to_stop=true;
	if(ite>0 && last<=sum_err1)to_stop=true;
	last=sum_err1;
	sum_err0=0;
	tot0=0;
	sum_err1=0;
	tot1=0;
}
}

void step(){
static int stepp;
cerr<<++stepp<<endl;
// cerr<<A.n<<endl;
// cerr<<rcnt<<endl;

	CUDA::radix_sort();
// exit(0);
	CUDA::forP<CUDA::upda1>();
	void*pp;
	cudaGetSymbolAddress(&pp,CUDA::has_contact);
	cudaMemset(pp,0,M*sizeof(bool));
	cudaGetSymbolAddress(&pp,CUDA::has_fluid);
	cudaMemset(pp,0,M*sizeof(bool));
	CUDA::forR<CUDA::upda2>();
	static const dim3 grid=dim3(1,CUDA::MY,CUDA::MZ),block=dim3(16,16);
	CUDA::for_neighbor<CUDA::update_rho,CUDA::update_rho_init,CUDA::update_rho_update,CUDA::update_rho_finish><<<grid,block>>>();
	CUDA::for_neighbor<CUDA::update2,CUDA::update2_init,CUDA::update2_update,CUDA::update2_finish><<<grid,block>>>();
	CUDA::for_neighbor<CUDA::update3,CUDA::update3_init,CUDA::update3_update,CUDA::update3_finish><<<grid,block>>>();
	CUDA::forR<CUDA::upda3>();
	CUDA::init_ite<<<1,1>>>();
	for(int ite=0;ite<10;++ite){
		CUDA::forR<CUDA::upda4>();
		bool flag_host;
		cudaMemcpyFromSymbol(&flag_host,CUDA::flag,1);
		if(!flag_host)CUDA::for_neighbor<CUDA::update_fluid_pressure_force,CUDA::update_fluid_pressure_force_init,CUDA::update_fluid_pressure_force_update,CUDA::update_fluid_pressure_force_finish><<<grid,block>>>();
		CUDA::forP<CUDA::upda5>();
		CUDA::for_neighbor<CUDA::update_s,CUDA::update_s_init,CUDA::update_s_update,CUDA::update_s_finish><<<grid,block>>>();
		CUDA::for_neighbor<CUDA::update_rigid_pressure_force,CUDA::update_rigid_pressure_force_init,CUDA::update_rigid_pressure_force_update,CUDA::update_rigid_pressure_force_finish><<<grid,block>>>();
		// cudaMemcpyFromSymbol(FR+1,CUDA::FR,rcnt*12);
		// cudaMemcpyFromSymbol(TR+1,CUDA::TR,rcnt*12);
		CUDA::forR<CUDA::upda6>();
		CUDA::forP<CUDA::upda7>();
		CUDA::forP<CUDA::upda8>();
		CUDA::for_neighbor<CUDA::update4,CUDA::update4_init,CUDA::update4_update,CUDA::update4_finish><<<grid,block>>>();
		CUDA::forP<CUDA::upda9>();
		CUDA::stop_ite<<<1,1>>>(ite);
		bool to_stop_host;
		cudaMemcpyFromSymbol(&to_stop_host,CUDA::to_stop,1);
		if(to_stop_host)break;
	}
	CUDA::forR<CUDA::upda4>();
	CUDA::for_neighbor<CUDA::update_fluid_pressure_force,CUDA::update_fluid_pressure_force_init,CUDA::update_fluid_pressure_force_update,CUDA::update_fluid_pressure_force_finish><<<grid,block>>>();
	CUDA::forP<CUDA::upda5>();
	CUDA::for_neighbor<CUDA::update_rigid_pressure_force,CUDA::update_rigid_pressure_force_init,CUDA::update_rigid_pressure_force_update,CUDA::update_rigid_pressure_force_finish><<<grid,block>>>();
	CUDA::forP<CUDA::upda7>();
	CUDA::forR<CUDA::upda10>();
	CUDA::forP<CUDA::upda11>();
	CUDA::forR<CUDA::upda12>();
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
	
	A.n=0;
	rcnt=0;
	A.move1();
	cudaMemcpyToSymbol(CUDA::A, &A, sizeof(Particles));
	cudaMemcpyToSymbol(CUDA::ri, ri+1, rcnt*sizeof(Rigid));
	cudaMemcpyToSymbol(CUDA::rcnt, &rcnt, sizeof(int));

	float tt;
	float rt=0;
	scanf("%d%f",&FC,&tt);
	for(int fr=0;fr<FC;){
		rt-=dt;
		if(rt<=0){
			cudaMemcpyFromSymbol(&A, CUDA::A, sizeof(Particles));
			A.move2();
			cudaMemcpyFromSymbol(ri+1,CUDA::ri,rcnt*sizeof(Rigid));
			
			for(;A.n<PC && time_p[A.n+1]<=fr;++A.n);
			for(;rcnt<RC && time_r[rcnt+1]<=fr;++rcnt);
			
			A.move1();
			cudaMemcpyToSymbol(CUDA::A, &A, sizeof(Particles));
			cudaMemcpyToSymbol(CUDA::ri, ri+1, rcnt*sizeof(Rigid));
			cudaMemcpyToSymbol(CUDA::rcnt, &rcnt, sizeof(int));
			A.move2();

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