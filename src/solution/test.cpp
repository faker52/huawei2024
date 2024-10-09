#include <windows.h>
#include <iostream>
#include <string>

using namespace std;

// 算法文件生成exe的路径
const string S_EXE = R"(..\cmake-build-debug\solution_main.exe)";
// 算法文件的路径
const string S_CPP = R"(..\src\solution\solution_main.cpp)";
// check_main.exe的路径
const string C_EXE = R"(..\cmake-build-debug\check_main.exe)";
// check_main.cpp的路径
const string C_CPP = R"(..\src\solution\check_main.cpp)";

/**
 * 从管道中读取数据，直到没有更多的数据为止
 */
string readPipeData(HANDLE hPipe) {
    const size_t bufferSize = 65536;  // 缓冲区大小，65536字节
    char buffer[bufferSize];
    string result;
    DWORD bytesRead;

    while (ReadFile(hPipe, buffer, bufferSize - 1, &bytesRead, NULL) || GetLastError() == ERROR_MORE_DATA) {
        buffer[bytesRead] = '\0'; // Ensure null termination
        result += buffer;
        if (bytesRead < bufferSize - 1) {
            break; // No more data
        }
    }

    return result;
}

/**
 * 编译
 */
void compileSource(const string& sourceFile, const string& outputFile) {
    string command = "g++ -o " + outputFile + " " + sourceFile;
    int result = system(command.c_str());
    if (result != 0) {
        cerr << "Compilation failed for " << sourceFile << "!" << endl;
        exit(1);
    }
}

int main() {
    const size_t bufferSize = 65536; // 缓冲区大小
    HANDLE hPipe1Read, hPipe1Write, hPipe2Read, hPipe2Write;
    PROCESS_INFORMATION pi1, pi2;
    STARTUPINFO si1, si2;
    SECURITY_ATTRIBUTES sa = { sizeof(SECURITY_ATTRIBUTES), NULL, TRUE };

    // 编译
    compileSource(C_CPP, C_EXE);
    compileSource(S_CPP, S_EXE);

    // 创建管道，hPipe1Read和hPipe1Write是第一个管道的读写句柄
    if (!CreatePipe(&hPipe1Read, &hPipe1Write, &sa, bufferSize)) {
        cerr << "CreatePipe failed" << endl;
        return 1;
    }

    // 创建管道，hPipe2Read和hPipe2Write是第二个管道的读写句柄。
    if (!CreatePipe(&hPipe2Read, &hPipe2Write, &sa, bufferSize)) {
        cerr << "CreatePipe failed" << endl;
        CloseHandle(hPipe1Read);
        CloseHandle(hPipe1Write);
        return 1;
    }

    // 启动 file1
    ZeroMemory(&si1, sizeof(si1));
    si1.cb = sizeof(si1);
    si1.hStdInput = hPipe2Read;
    si1.hStdOutput = hPipe1Write;
    si1.hStdError = GetStdHandle(STD_ERROR_HANDLE);
    si1.dwFlags |= STARTF_USESTDHANDLES;
    ZeroMemory(&pi1, sizeof(pi1));
    if (!CreateProcess(NULL, const_cast<char*>(S_EXE.c_str()), NULL, NULL, TRUE, 0, NULL, NULL, &si1, &pi1)) {
        cerr << "CreateProcess for file1 failed" << endl;
        CloseHandle(hPipe1Read);
        CloseHandle(hPipe1Write);
        CloseHandle(hPipe2Read);
        CloseHandle(hPipe2Write);
        return 1;
    }

    // 启动 file2
    ZeroMemory(&si2, sizeof(si2));
    si2.cb = sizeof(si2);
    si2.hStdInput = hPipe1Read;
    si2.hStdOutput = hPipe2Write;
    si2.hStdError = GetStdHandle(STD_ERROR_HANDLE);
    si2.dwFlags |= STARTF_USESTDHANDLES;
    ZeroMemory(&pi2, sizeof(pi2));
    if (!CreateProcess(NULL, const_cast<char*>(C_EXE.c_str()), NULL, NULL, TRUE, 0, NULL, NULL, &si2, &pi2)) {
        cerr << "CreateProcess for file2 failed" << endl;
        TerminateProcess(pi1.hProcess, 1);
        CloseHandle(hPipe1Read);
        CloseHandle(hPipe1Write);
        CloseHandle(hPipe2Read);
        CloseHandle(hPipe2Write);
        return 1;
    }

    // 关闭不需要的句柄
    CloseHandle(hPipe1Write);
    CloseHandle(hPipe2Write);

    // 等待进程结束
    WaitForSingleObject(pi1.hProcess, INFINITE);
    WaitForSingleObject(pi2.hProcess, INFINITE);

    // 获取 file1 和 file2 的退出代码
    DWORD exitCodeFile1, exitCodeFile2;
    if (!GetExitCodeProcess(pi1.hProcess, &exitCodeFile1)) {
        cerr << "Failed to get exit code for file1" << endl;
        return 1;
    }
    if (!GetExitCodeProcess(pi2.hProcess, &exitCodeFile2)) {
        cerr << "Failed to get exit code for file2" << endl;
        return 1;
    }

    // 输出退出代码
    cout << S_EXE + " exit code: " << exitCodeFile1 << endl;
    cout << C_EXE + " exit code: " << exitCodeFile2 << endl;

    // 读取 solution_main.exe 的输出
    string solutionMainOutput = readPipeData(hPipe1Read);
    cout << S_EXE +  " output:\n" << solutionMainOutput << endl;

    // 读取 check_main.exe 的输出
    string checkMainOutput = readPipeData(hPipe2Read);
    cout << C_EXE + " output:\n" << checkMainOutput << endl;

    // 清理句柄
    CloseHandle(hPipe1Read);
    CloseHandle(hPipe2Read);
    CloseHandle(pi1.hProcess);
    CloseHandle(pi1.hThread);
    CloseHandle(pi2.hProcess);
    CloseHandle(pi2.hThread);

    return 0;
}
