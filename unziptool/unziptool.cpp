// unziptool.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include<Windows.h>
#include<tchar.h>
#include<atlconv.h>
#include<string>
#include <iostream>
#include<string>
#include<algorithm>
#include"zip.h"
#include"unzip.h"
using namespace std;
int main(int argc, char* argv[])
{
    if (argc != 3) return -1;
    string srcFile(argv[1]);
    string dstPath(argv[2]);
    USES_CONVERSION;
    std::replace(srcFile.begin(), srcFile.end(), '/', '\\');
    HZIP hz = OpenZip(A2W(srcFile.c_str()), 0);
    ZIPENTRY ze; GetZipItem(hz, -1, &ze); int numitems = ze.index;
    for (int i = 0; i < numitems; i++)
    {
        GetZipItem(hz, i, &ze);
        std::string str = dstPath + string("\\") + std::string(W2A(ze.name));
        std::replace(str.begin(), str.end(), '/', '\\');
        UnzipItem(hz, i, A2W(str.c_str()));
    }
    CloseZip(hz);
    return 0;
}

// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
