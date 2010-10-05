# Microsoft Developer Studio Project File - Name="guess_dll_vc6" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Dynamic-Link Library" 0x0102

CFG=guess_dll_vc6 - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "guess_dll_vc6.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "guess_dll_vc6.mak" CFG="guess_dll_vc6 - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "guess_dll_vc6 - Win32 Release" (based on "Win32 (x86) Dynamic-Link Library")
!MESSAGE "guess_dll_vc6 - Win32 Debug" (based on "Win32 (x86) Dynamic-Link Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
MTL=midl.exe
RSC=rc.exe

!IF  "$(CFG)" == "guess_dll_vc6 - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "ReleaseDLL"
# PROP Intermediate_Dir "ReleaseDLL"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /MT /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "GUESS_DLL_VC6_EXPORTS" /YX /FD /c
# ADD CPP /nologo /MT /W3 /GX /O2 /I "../framework" /I "../modules" /I "../libraries/gutil" /I "../libraries/plib" /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "GUESS_DLL_VC6_EXPORTS" /YX /FD /c
# ADD BASE MTL /nologo /D "NDEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "NDEBUG" /mktyplib203 /win32
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /dll /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /dll /machine:I386 /out:"ReleaseDLL/guess.dll"

!ELSEIF  "$(CFG)" == "guess_dll_vc6 - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "DebugDLL"
# PROP Intermediate_Dir "DebugDLL"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /MTd /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "GUESS_DLL_VC6_EXPORTS" /YX /FD /GZ /c
# ADD CPP /nologo /MTd /W3 /Gm /GX /ZI /Od /I "../framework" /I "../modules" /I "../libraries/gutil" /I "../libraries/plib" /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "GUESS_DLL_VC6_EXPORTS" /YX /FD /GZ /c
# ADD BASE MTL /nologo /D "_DEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "_DEBUG" /mktyplib203 /win32
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /dll /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /dll /debug /machine:I386 /out:"DebugDLL/guess.dll" /pdbtype:sept

!ENDIF 

# Begin Target

# Name "guess_dll_vc6 - Win32 Release"
# Name "guess_dll_vc6 - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\modules\canexch.cpp
# End Source File
# Begin Source File

SOURCE=..\modules\driver.cpp
# End Source File
# Begin Source File

SOURCE=..\modules\growth.cpp
# End Source File
# Begin Source File

SOURCE=..\framework\guess.cpp
# End Source File
# Begin Source File

SOURCE=..\modules\guessio.cpp
# End Source File
# Begin Source File

SOURCE=..\cru\guessio\guessio_cru.cpp
# End Source File
# Begin Source File

SOURCE=..\libraries\gutil\gutil.cpp
# End Source File
# Begin Source File

SOURCE=..\windows_version\main.cpp
# End Source File
# Begin Source File

SOURCE=..\libraries\plib\plib.cpp
# End Source File
# Begin Source File

SOURCE=..\modules\soilwater.cpp
# End Source File
# Begin Source File

SOURCE=..\modules\somdynam.cpp
# End Source File
# Begin Source File

SOURCE=..\modules\vegdynam.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\modules\canexch.h
# End Source File
# Begin Source File

SOURCE=..\framework\config.h
# End Source File
# Begin Source File

SOURCE=..\cru\guessio\cru_1901_2006.h
# End Source File
# Begin Source File

SOURCE=..\cru\guessio\cru_1901_2006misc.h
# End Source File
# Begin Source File

SOURCE=..\modules\driver.h
# End Source File
# Begin Source File

SOURCE=..\modules\growth.h
# End Source File
# Begin Source File

SOURCE=..\framework\guess.h
# End Source File
# Begin Source File

SOURCE=..\modules\guessio.h
# End Source File
# Begin Source File

SOURCE=..\libraries\gutil\gutil.h
# End Source File
# Begin Source File

SOURCE=..\windows_version\main.h
# End Source File
# Begin Source File

SOURCE=..\libraries\plib\plib.h
# End Source File
# Begin Source File

SOURCE=..\modules\soilwater.h
# End Source File
# Begin Source File

SOURCE=..\modules\somdynam.h
# End Source File
# Begin Source File

SOURCE=..\modules\vegdynam.h
# End Source File
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# End Target
# End Project
