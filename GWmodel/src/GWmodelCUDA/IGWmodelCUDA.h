#pragma once

#ifdef WIN32
#ifdef CREATDLL_EXPORTS
#define GWMODELCUDA_API __declspec(dllexport)
#else
#define GWMODELCUDA_API __declspec(dllimport)
#endif // CREATDLL_EXPORTS
#else
#define GWMODELCUDA_API 
#endif

class GWMODELCUDA_API IGWmodelCUDA 
// 刚刚这里好像写反了，不知道是不是这个原因
// 好像是可以了
// 这里咋会反过来了 笑死
// 还有 Makevars.win 要不要给 32 位的加一下，xm
// 加了的话，就得再弄dll，文件夹会巨大无比
// 我们把这个 DLL 和 LIB 放到 Github 上，不放到包里。就是包里只提供源码，如果是Windows下想用 CUDA,就自己编译，我们负责提供 DLL, lib
// 你说起这个倒是有点意思，我看看Rtools里面有
// 你是说 CURL？
// 嗯对
// 可以，很强，如果有这个就好了，放在 github ，编译的时候自动下载。但是 Windows 能用 configure 么
// 没移植， 那个windows下有没有啥命令是类似于curl，wget的？windows自带的， 可以用configure，但是必须要移植了的命令才可以写
// powershell 倒是有，CMD 没有
// 哪个命令？
// $client = new-object System.Net.WebClient
// $client.DownloadFile('URL', '保存位置/文件名称.文件类型')
// 这样可以下载文件
// 这是VB？主要是试试能不能用Rtools里面的sh命令执行他
// .NET 系列的，powershell 的语法，里面可以用很多 .NET 的对象
// 
{
public:
	virtual void SetX(int i, int k, double value) = 0;
	virtual void SetY(int i, double value) = 0;
	virtual void SetDp(int i, double u, double v) = 0;
	virtual void SetRp(int i, double u, double v) = 0;
	virtual void SetDmat(int i, int j, double value) = 0;

	virtual double GetBetas(int i, int k) = 0;
	virtual double GetBetasSE(int i, int k) = 0;
	virtual double GetShat1() = 0;
	virtual double GetShat2() = 0;
	virtual double GetQdiag(int i) = 0;


	virtual bool Regression(
		bool hatmatrix,
		double p, double theta, bool longlat,
		double bw, int kernel, bool adaptive,
		int groupl, int gpuID
	) = 0;

	virtual double CV(
		double p, double theta, bool longlat,
		double bw, int kernel, bool adaptive,
		int groupl, int gpuID
	) = 0;
};

extern "C" GWMODELCUDA_API IGWmodelCUDA* GWCUDA_Create(int N, int K, bool rp_given, int n, bool dm_given);
extern "C" GWMODELCUDA_API void GWCUDA_Del(IGWmodelCUDA* pInstance);