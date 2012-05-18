////////////////////////////////////////////////////
// Plugin.h
// Standard header file for all AmiBroker plug-ins
//
// Copyright (C)2001 Tomasz Janeczko, amibroker.com
// All rights reserved.
// 
// Last modified: 2001-09-24 TJ
// 
// You may use this code in your own projects provided that:
//
// 1. You are registered user of AmiBroker
// 2. The software you write using it is for personal, noncommercial use only
//
// For commercial use you have to obtain a separate license from Amibroker.com
//
////////////////////////////////////////////////////

#ifndef PLUGIN_H
#define PLUGIN_H 1

// In MSVC you must add -D_USRDLL when building the DLL, and not when using it.
// This is because all exportable functions must have undecorated names 
// Under Unix/Linux, this is automatically ignored.
#if    defined(__linux) || defined(__unix)
#define PLUGINAPI
#else
#ifdef _USRDLL
#define PLUGINAPI extern "C" __declspec(dllexport)
#else
#define PLUGINAPI extern "C" __declspec(dllimport)
#endif
#endif

// define signed and unsigned one byte types
typedef unsigned char UBYTE;
typedef signed char SBYTE;

// useful macros for empty values
#define EMPTY_VAL (-1e10f)
#define IS_EMPTY( x ) ( x == EMPTY_VAL )
#define NOT_EMPTY( x ) ( x != EMPTY_VAL )

// the list of AmiVar types
enum { VAR_NONE, VAR_FLOAT, VAR_ARRAY, VAR_STRING, VAR_DISPATCH };

// AmiVar is a variant-like structure/union
// that holds any AFL value
// type member holds variable type (see VAR_ enum above)

typedef struct AmiVar
{
    int type;
    union 
    {
        float   val;
        float   *array;
        char    *string;
		void	*disp;
    };
} AmiVar;


// PluginInfo structure holds
// general information about plugin 
struct PluginInfo
{
	int			nStructSize;	// this is sizeof( struct PluginInfo )
	int			nType;			// plug-in type currently 1
								// [indicator is the only one supported]
	int			nVersion;		// plug-in version coded to int as
								// MAJOR*10000 + MINOR * 100 + RELEASE
	int			nFlags;			// for future use - set to 0
	char		szName[ 64 ];	// long name displayed in the Plugin dialog
	char    	szVendor[ 64 ];	// name of the plug-in vendor
	int			nCertificate;	// certificate code - zero for private plug-ins
	int			nMinAmiVersion;	// minimum required AmiBroker version
								// (should be >= 380000 -> AmiBroker 3.8)
};

// SiteInterface structure
// Defines call-back function pointers. The structure is filled with correct
// pointers by the AmiBroker and passed to DLL via SetSiteInterface() call.
// 
// SiteInterface is used as a way to call-back AmiBroker built-in functions
//

struct SiteInterface 
{
	int			nStructSize;
	int			(*GetArraySize) (void);	   
	float*		(*GetStockArray)( int nType );
	AmiVar		(*GetVariable) ( const char *pszName );
	void		(*SetVariable) ( const char *pszName, AmiVar newValue );
	AmiVar		(*CallFunction) ( const char *szName,
								  int nNumArgs,
								  AmiVar *ArgsTable );
	AmiVar		(*AllocArrayResult) (void);
	void *		(*Alloc) (unsigned int nSize);
	void		(*Free) (void *pMemory);
};

// FunDesc structure
// Holds the pointer to user-defined function that can be called by AmiBroker.
// It holds also the number of array, string, float and default arguments
// for the function and the default values.
//
typedef struct FunDesc
{
    AmiVar (*Function)( int NumArgs, AmiVar *ArgsTable );
    UBYTE   ArrayQty;       // number of Array arguments required   
    UBYTE   StringQty;      // number of String arguments required
    SBYTE   FloatQty;       // number of float args 
    UBYTE   DefaultQty;     // number of default float args
    float   *DefaultValues; // the pointer to defaults table 
} FunDesc;


// FunctionTag struct
// Holds the Name of the function and the corresponding FunDesc structure.
// This structure is used to define function table that is retrieved by
// AmiBroker via GetFunctionTable() call when AFL engine is initialized.
// That way new function names are added to the AFL symbol table
// and they become accessible.

typedef struct FunctionTag
{
	char	*Name;
	FunDesc	 Descript;
} FunctionTag;


///////////////////////////////////////////////////
// EXPORTED FUNCTONS
//
// Each AmiBroker plug-in DLL must export the following
// functions:
// 1. GetPluginInfo	- called when DLL is loaded
// 2. Init - called when AFL engine is being initialized 
// 3. Release - called when AFL engine is being closed
// 4. GetFunctionTable - called when AFL engine is being initialized 
// 5. SetSiteInteface - called when AFL engine is being initialized 
//
// Each function may be called multiple times.
//
// The order of calling functions during intialization is
// as follows:
//
// SetSiteInterface -> GetFunctionTable	-> Init -> 
// ... normal work ....
// Release
//
// This cycle may repeat multiple times
// 
// All functions in the plug in DLL use _cdecl calling convention
// (the default for C compiler)

PLUGINAPI int GetPluginInfo( struct PluginInfo *pInfo );

PLUGINAPI int Init(void);
PLUGINAPI int Release(void);

PLUGINAPI int GetFunctionTable( FunctionTag **ppFunctionTable );
PLUGINAPI int SetSiteInterface( struct SiteInterface *pInterface );

////////////////////////////////////////
// Global-scope data
////////////////////////////////////////

// FunctionTable should be defined 
// in the implementation file of your functions
extern FunctionTag gFunctionTable[];
extern int		   gFunctionTableSize;

// Site interface is defined in Plugin.cpp
extern struct SiteInterface gSite;


#endif
