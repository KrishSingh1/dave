///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                  PLIB Version 2.1                                 //
//                              (Fully portable version)                             //
//                                                                                   //
//                                Written by Ben Smith                               //
//                                 University of Lund                                //
//                                                                                   //
//  PLIB is a utility library for reading input from a profile script. This is a     //
//  text file consisting of commands and data, in a format defined in part by PLIB,  //
//  and in part by the calling program. Language elements supported are:             //
//                                                                                   //
//  1. COMMANDS: statements beginning with an identifier (a string of non-white      //
//     space characters commencing with an alphabetic character and not including    //
//     any of the following: - + . " ' ( ) !. The identifier may be followed by      //
//     a string (enclosed in single or double quotation marks) with a specified      //
//     maximum significant length, a specified number of integers or real numbers    //
//     in a specified range, a specified number of bool values (numerals             //
//     interpreted such that 0=false, non-zero=true) or no parameters at all.        //
//     Identifiers are non-case-sensitive.                                           //
//     e.g. title "Kiruna"                                                           //
//          mtemp -14.2 -13.5 -9.5 -3.5 2.9 9.3 12.6 10.3 4.8 -1.5 -8.1 -12.1        //
//          include 1                                                                //
//          exit                                                                     //
//                                                                                   //
//  2. SETS: named blocks of script enclosed in parentheses. Sets may be nested.     //
//     e.g. taxon "picea" ( ... )                                                    //
//                                                                                   //
//  3. GROUPS: named blocks of script enclosed in parentheses which are interpreted  //
//     as "macros" when the script is interpreted. The contents of the group are     //
//     implicitly inserted in the script wherever the name of the group appears      //
//     following the group definition. Groups may be nested, but always have         //
//     'global' scope. Groups are identified by the keyword group followed by the    //
//     group name in single or double quotation marks.                               //
//     e.g. group 'woody' ( tree 1 wooddens 250 leaftoroot 1.0 )                     //
//          taxon "picea" ( woody )                                                  //
//                                                                                   //
//  4. COMMENTS: strings of characters beginning with an exclamation mark (!). They  //
//     are ignored to the end of the line they appear on.                            //
//     e.g. !this is a comment                                                       //
//                                                                                   //
//  -------------------------------------------------------------------------------  //
//                         LINKING TO PLIB, AND DEPENDENCIES                         //
//  -------------------------------------------------------------------------------  //
//                                                                                   //
//  The recommended way of building an executable incorporating the functionality    //
//  of PLIB is to first build PLIB as a binary library and link this library to      //
//  the calling program. Each source code file containing calls to PLIB library      //
//  functions must refer (via a suitable #include directive) to the PLIB header      //
//  file plib.h.                                                                     //
//                                                                                   //
//  Character strings are represented in PLIB as xtring objects. String arguments    //
//  to PLIB library functions are also of type xtring. Therefore, both PLIB and the  //
//  calling program must be linked to the Xtring library when built. Further         //
//  information regarding the Xtring library is available from the author (see       //
//  below).                                                                          //
//                                                                                   //
//  -------------------------------------------------------------------------------  //
//                      READING INFORMATION FROM A PLIB SCRIPT                       //
//  -------------------------------------------------------------------------------  //
//                                                                                   //
//  To process a PLIB script, call function plib, specifying the full pathname of    //
//  the script file to process:                                                      //
//                                                                                   //
//  bool plib(xtring filename)                                                       //
//                                                                                   //
//  Function plib returns TRUE if there were no errors, otherwise FALSE (zero).      //
//  In the event of an error, a explanatory message is sent to function              //
//  plib_receivemessage in the calling program (see below).                          //
//                                                                                   //
//  -------------------------------------------------------------------------------  //
//                COMMUNICATION BETWEEN PLIB AND THE CALLING PROGRAM                 //
//  -------------------------------------------------------------------------------  //
//                                                                                   //
//  Communication is via PLIB library functions, accessible to the calling program.  //
//  The calling program must also implement two, and in some cases, three, special   //
//  functions accessible to PLIB.                                                    //
//                                                                                   //
//  PUBLIC FUNCTIONS IN THE PLIB LIBRARY (see below):                                //
//     bool plib(xtring filename)                                                    //
//     declareitem([various overloads])                                              //
//     void callwhendone(int callback)                                               //
//     bool itemparsed(xtring identifier)                                            //
//     void sendmessage(xtring header,xtring message)                                //
//     void plibhelp()                                                               //
//     void plibabort()                                                              //
//                                                                                   //
//  FUNCTIONS TO BE IMPLEMENTED BY THE CALLING PROGRAM (see below):                  //
//     void plib_declarations(int id,xtring setname)                                 //
//     void plib_receivemessage(xtring text)                                         //
//     void plib_callback(int callback)                                              //
//                                                                                   //
//  The calling program must supply functions plib_declarations and                  //
//  plib_receivemessage and, if the callback feature is implemented for any command  //
//  (see below), function plib_callback.                                             //
//                                                                                   //
//  Function plib_declarations is called by PLIB just before interpretation of a     //
//  script begins, and subsequently whenever a set header is encountered. It should  //
//  include calls to the functions declareitem which define the identifiers and      //
//  formats for commands and sets in the script (or this set), and what action       //
//  should be taken when a particular command or set identifier is encoutered. The   //
//  following alternative forms of function declareitem are available:               //
//                                                                                   //
//  declareitem(xtring identifier,xtring* param,int maxlen,int callback)             //
//    When the string specified by identifier is encountered, expect the next item   //
//    to be a string and write it to 'param'. 'maxlen' specifies the maximum         //
//    allowable length of the string, not including the terminating null character;  //
//    the string is truncated, if necessary, before being written to 'param'. If     //
//    'callback' is non-zero, function plib_callback (see below) is called with the  //
//    integer value of 'callback' as a parameter, AFTER the assignment to param.     //
//    e.g. title "kiruna"                                                            //
//                                                                                   //
//  declareitem(xtring identifier,int* param,int min,int max,int nparam,             //
//              int callback)                                                        //
//    Here the number of integers specified in 'nparam' are expected after the       //
//    identifier and assigned to the array (or simple variable) pointed to by        //
//    'param'. Floating-point numbers are rounded to the nearest integer. The        //
//    parameters should be in the range 'min'-'max', otherwise a PLIB error          //
//    results.                                                                       //
//    e.g. npat 1000                                                                 //
//                                                                                   //
//  declareitem(xtring identifier,double* param,double min,double max,int nparam,    //
//              int callback)                                                        //
//    Here the number of floating-point numbers specified in 'nparam' are expected   //
//    after the identifier and assigned to the array (or simple variable) pointed    //
//    to by 'param'. The parameters should be in the range 'min'-'max', otherwise    //
//    a PLIB error results.                                                          //
//    e.g. mtemp -14.2 -13.5 -9.5 -3.5 2.9 9.3 12.6 10.3 4.8 -1.5 -8.1 -12.1         //
//                                                                                   //
//  declareitem(xtring identifier,bool* param,int nparam,int callback)               //
//    Here the number of numbers (nominally integers 0 or 1) specified in 'nparam'   //
//    are expected after the identifier and are reinterpreted as bools (0=false,     //
//    non-zero=true) and assigned to the array (or simple variable) pointed to by    //
//    'param'.                                                                       //
//    e.g. include 1                                                                 //
//                                                                                   //
//  declareitem(xtring identifier,bool* param,int callback)                          //
//    If this identifier is encountered, the variable pointed to by 'param' is       //
//    assigned the value true.                                                       //
//    e.g. exit                                                                      //
//                                                                                   //
//  declareitem(xtring identifier,int id,int callback)                               //
//    Specifies the identifier and associated id-code for a set header. If the       //
//    specified identifier is encountered, it is interpreted as a set header and     //
//    should be followed by the set name as a string and an open parenthesis.        //
//    Function plib_declarations is then called, passing 'id' from this call to      //
//    declareitem and the set name encountered in the script.                        //
//    e.g. taxon "picea" (                                                           //
//                                                                                   //
//  The following additional overloads of declareitem are available, and include     //
//  argument 'xtring help,' a string of text up to 62 characters in length           //
//  (characters in excess of this are ignored) which is output if function plibhelp  //
//  is called (further information below). In other respects the following forms     //
//  have identical function to the corresponding forms excluding argument            //
//  'xtring help', and listed above.                                                 //
//                                                                                   //
//  declareitem(xtring identifier,xtring* param,int maxlen,int callback,             //
//              xtring help)                                                         //
//  declareitem(xtring identifier,int* param,int min,int max,int nparam,             //
//              int callback,xtring help)                                            //
//  declareitem(xtring identifier,double* param,double min,double max,int nparam,    //
//              int callback,xtring help)                                            //
//  declareitem(xtring identifier,bool* param,int nparam,int callback,               //
//              xtring help)                                                         //
//  declareitem(xtring identifier,bool* param,int callback,xtring help)              //
//  declareitem(xtring identifier,int id,int callback,xtring help)                   //
//                                                                                   //
//  Functions declareitem are of type bool and return false if dynamic memory        //
//  could not be allocated, otherwise true.                                          //
//                                                                                   //
//  Function callwhendone declares a callback code which, if non-zero, is passed in  //
//  a call to function plib_callback following input of the last recognised PLIB     //
//  item in the current set. This can provide the calling program with an            //
//  opportunity to query, via function itemparsed, which of the items declared for   //
//  this set were actually encountered in the script.                                //
//                                                                                   //
//  Function itemparsed should be called from function plib_callback following       //
//  input of data for a particular set (see function callwhendone above). It         //
//  returns true if the string 'identifier' corresponds to one of the PLIB           //
//  commands declared by calls to function declareitem for this set, AND this        //
//  command was encountered in processing this set in the script. If 'identifier'    //
//  is not recognised, or was not encountered in the script, the return value is     //
//  false.                                                                           //
//                                                                                   //
//  Function plib_receivemessage is called by PLIB whenever output is required to    //
//  be conveyed to the user. Plib sends output in the form of an xtring (character   //
//  string) object. It is the responsibility of the calling program to send this     //
//  string to an output stream or device (typically the screen, or a text file).     //
//  Output generated by PLIB includes error messages, help text and messages         //
//  originating from the calling program and sent to PLIB via function sendmessage.  //
//                                                                                   //
//  Function sendmessage may be called by the calling program to output warning or   //
//  informational messages to the user during processing of a PLIB script.           //
//  Typically sendmessage would be called from within function plib_callback.        //
//  Parameter 'heading' should consist of a keyword identifying the type of message  //
//  sent (e.g. "Warning"). More detailed information may be given within parameter   //
//  'message'.                                                                       //
//                                                                                   //
//  Documentation of a calling module's native keywords (as defined in function      //
//  plib_declarations) is possible using the 'xtring help' overloads of function     //
//  declareitem (above). A list of all keywords and their associated help text (if   //
//  provided) is sent as sequential calls (one per line of output text) to function  //
//  plib_receivemessage, when function plibhelp is called from the calling program.  //
//                                                                                   //
//  Finally, function plibabort is provided to allow processing of a PLIB script to  //
//  be forcibly aborted. Typically, plibabort would be called from within function   //
//  plib_callback.                                                                   //
//                                                                                   //
//  -------------------------------------------------------------------------------  //
//  This version dated 5 November 2004                                               //
//  E-mail: ben@planteco.lu.se                                                       //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

class xtring;

// PLIB LIBRARY FUNCTIONS

bool plib(xtring filename);
bool declareitem(xtring identifier,xtring* param,int maxlen,int callback);
bool declareitem(xtring identifier,int* param,int min,int max,int nparam,
	int callback);
bool declareitem(xtring identifier,double* param,double min,double max,int nparam,
	int callback);
bool declareitem(xtring identifier,bool* param,int nparam,int callback);
bool declareitem(xtring identifier,bool* param,int callback);
bool declareitem(xtring identifier,int id,int callback);
bool declareitem(xtring identifier,xtring* param,int maxlen,int callback,
	xtring help);
bool declareitem(xtring identifier,int* param,int min,int max,int nparam,
	int callback,xtring help);
bool declareitem(xtring identifier,double* param,double min,double max,int nparam,
	int callback,xtring help);
bool declareitem(xtring identifier,bool* param,int nparam,int callback,
	xtring help);
bool declareitem(xtring identifier,bool* param,int callback,xtring help);
bool declareitem(xtring identifier,int id,int callback,xtring help);
void callwhendone(int callback);
bool itemparsed(xtring identifier);
void sendmessage(xtring heading,xtring message);
void plibhelp();
void plibabort();

// FUNCTIONS TO BE IMPLEMENTED BY THE CALLING PROGRAM

void plib_declarations(int id,xtring setname);
void plib_callback(int callback);
void plib_receivemessage(xtring text);
