///////////////////////////////////////////////////////////////////////////////////////
// MODULE SOURCE CODE FILE
//
// Module:                Main module for Linus/Unix command line version of LPJ-GUESS
//                        Includes support for X11 graphical output
//                        Output from global functions plot(), resetwindow() and
//                        clear_all_graphs() are sent to an X Windows server if
//                        available
// Header file name:      main.h
// Source code file name: main.cpp
// Written by:            Ben Smith
// Version dated:         2004-08-10

#include "config.h"
#include "guess.h"
#include <stdarg.h>
#include <X11/Xlib.h>
#include <X11/Xatom.h>
#include <unistd.h>
#include <pthread.h>


///////////////////////////////////////////////////////////////////////////////////////
// GLOBAL SETTINGS
// Some can be changed from the default by options at command line

char WSFILE[]="winspex.bin";
	// Name of file to save window dimensions between runs
bool iflog=true;
	// Whether to output to a log file
	// Change by command line option -nolog
xtring file_log="guess.log";
	// Name of log file to which output from all dprintf and fail calls is sent
	// Change by command line option -log <filename>
bool ifgra=true;
	// Whether to perform graphical output
	// Change by command line option -nogra
bool ifwait=false;
	// Whether to wait for Enter from stdin at end of run
	// (to preserve graphics)
bool ifmono=false;
	// Whether graphical output should be monochrome (default colour)
	// Change by command line option -mono
xtring display="";
	// Address/ip-number to X Windows manager
	// null implies DISPLAY environmental variable used
	// Override by command line option -display <address>
int width=320,height=240;
	// Width and height of graphics window in pixels
	// Override by command line options -width <width> and -height <height>
bool use_saved_window_dims=true;
	// Whether saved window dimensions and position from a previous run
	// should override default width and height
	// Specifying -width or -height on the command line sets this to false

	

///////////////////////////////////////////////////////////////////////////////////////
// DECLARATION OF GRAPHICS MANAGER CLASS

class Gp; // (forward declaration)
class XWindow; // (forward declaration)
class XClient; // (forward declaration)

const double MAXBGAUTOTIC=6.0; // maximum number of tics on automatic scaling

// Line patterns for monochrome graphs

const double SCALE_LINE_SYMBOL=2.0;
	// Scaling factor for monochrome line symbol

const int NLINESYMBOL=10;
const int MAXSECTION=6;

struct Linesymbol {
	int nsection;
	int section[MAXSECTION];
};

const Linesymbol linesymbol[NLINESYMBOL]={
	1,{100},
	2,{1,1},
	2,{6,2},
	4,{6,2,1,2},
	6,{6,2,1,2,1,2},
	2,{2,6},
	2,{6,6},
	2,{10,2},
	4,{6,6,1,6},
	6,{6,6,1,2,1,6}
};


// Text justification specifiers
enum justtype{AUTOJUST,TL,TC,TR,CL,CC,CR,BL,BC,BR,LEFT,RIGHT,CENTRE};

struct Location {
	int h;
	int v;
};

struct Ticdef {

	// Tic format

	bool centre;
	int length;
};

class Axisdef {

	// Axis format

	friend class Gp;
	friend class CGraph;

	double displ;
	bool reverse; // if true, axis goes from top to bottom instead
	bool log; // if true, log transformed
	bool zero; // if true, axis starts at zero if possible

	Axisdef() {
		displ=0.0;
		reverse=false;
		log=false;
		zero=false;
	}
};

class Labeldef {
	
	// Format for tic labels

	friend class Gp;
	friend class CGraph;

	int width;
	int places;
	bool exp;

	Labeldef() {
		width=0;
		places=0;
		exp=false;
	}
};

struct Data {
	
	// Struct containing data value as part of a linked list

	Data* pnext;

	bool plotted;
	float x;
	float y;
};

class Dataseries {
	
	// Struct containing data series for a graph

	friend class Gp;
	friend class CGraph;

	int series; // id code for series
	char name[80]; // name of series
	Dataseries* pnext; // pointer to next series
	Data* pfirst; // pointer to first data item in series
	int colourid; // colour id number (from palette)

	Dataseries(int id,int colour,char* n) {
		series=id;
		pfirst=NULL;
		pnext=NULL;
		colourid=colour;
		strcpy(name,n);
	}
};

class Gp {
	
	// Class defining parameters for a particular graph (window)

	friend class CGraph;

public:
	Gp* pnext;
	char name[80];
	int tag;
	XWindow* pwnd;
	bool m_update;

	// Window width and height
	int width;
	int height;
	int margin;

	// Vertical and horizontal data area parameters
	int hmin;
	int hrange;
	int vmax;
	int vrange;
	int dir;

	// Scaling parameters for data area
	double xbase;
	double xmul;
	double ybase;
	double ymul;
	double xmin;
	double xmax;
	double ymin;
	double ymax;

	double scale;
	int textheight;

	// Major and minor tic intervals
	double autoticx;
	double minorticx;
	double autoticy;
	double minorticy;

	// Axis styles
	Axisdef xaxispars;
	Axisdef yaxispars;

	// Tic styles
	Ticdef ticstyle;
	Ticdef minorticstyle;
	bool wantminortics;

	// Axis label styles
	Labeldef xlabelstyle;
	Labeldef ylabelstyle;

	// amount of data entered so far
	int ndata;

	// Data series
	int nseries; 
	Dataseries* pfirst;

	// Legend parameters
	int leglinelen; // legend line length in pixels

public:
	void configure(int w,int h,int textsize) {
		width=w;
		height=h;
		hmin=width/7;
		hrange=width-width/7-width/10;
		vrange=height-height/10-height/15;
		vmax=height-height/10;
		ticstyle.length=max(4,min(w,h)*0.015);
		minorticstyle.length=ticstyle.length*0.5+1;
		leglinelen=ticstyle.length*8; // previously 5 BLARP!!
		scale=min(w,h)*0.00417;
		textheight=textsize;
	};

private:
	Gp(char* str,int w,int h,int textsize) {
		tag=0; // mark as untagged
		strcpy(name,str);
		pfirst=NULL;
		ticstyle.centre=false;
		minorticstyle.centre=false;
		wantminortics=true;
		ndata=0;
		nseries=0;
		dir=1; // vertical "direction"
		margin=0;
		configure(w,h,textsize);
	}

	// Translation of data values into window pixel positions and vice-versa
	int gtopx(double g) {
		return (int)((g-xbase)*xmul*(double)hrange+0.5)+hmin;
	}
	int gtopy(double g) {
		return vmax-(int)((g-ybase)*ymul*(double)vrange+0.5);
	}
	double ptogx(int g) {
		return (double)(g-hmin)/(double)hrange/xmul+xbase;
	}
	double ptogy(int g) {
		return ybase+((double)(vmax-g)/(double)vrange)/ymul;
	}
};


class CGraph {

	// Graphics manager class
	// Applications should maintain one CGraph object for each XClient and
	// X server connection

private:
	int nwindow;
	Gp* firstwindow;

	// MEMBER FUNCTIONS
private:
	void autoscaleaxis(Axisdef& axispars,Labeldef& labelstyle,double& base,double& mul,
		double& setmin,double& setmax,double& autotic,double& minortic,
		double min,double max,double round,bool fixrange);
	void setlabelstyle(Labeldef& labelstyle,int width,int places,bool exp);
	void bgsprint(char* string,double var,int width,int places,bool exp);
	void logprint(char* labtext,double ctic,bool exp);
	void pline(XWindow* pwin,int hstart,int vstart,int hend,int vend);
	void ptext(XWindow* pwin,Gp& gp,int hpos,int vpos,justtype just,char* text);
	void plotdata(XWindow* pwin,Gp& gp,Dataseries& series);
	void xaxis(XWindow* pwin,Gp& gp);
	void yaxis(XWindow* pwin,Gp& gp);
	void xtics(XWindow* pwin,Gp& gp);
	void ytics(XWindow* pwin,Gp& gp);
	void legend(XWindow* pwin,Gp& gp);
	Gp& getgp(char* name);

public:
	void cleanup();
	void cleargp(Gp* pgpx);
	void reset();
	void resetgp(Gp* pgp);
	void resetwindow(char* wname);

public:
	CGraph() {
		firstwindow=NULL;
		nwindow=0;
	}
	bool createwindow(XClient* pxclient,char* wname,unsigned int width,
		unsigned int height,bool use_saved_dimensions);
	void tag(int n);
	void scalex(char* name,double min,double max,double round);
	void scaley(char* name,double min,double max,double round);
	bool newseries(char* wname,char* name);
	bool data(char* wname,char* name,double xval,double yval,bool rescale);
	bool data(char* wname,char* name,double xval,double yval) {
		return data(wname,name,xval,yval,true);
	}
	void update(Gp* pgp);
	void save(FILE*& out,char* name);
	void destroyall() {
		Data* d,*dnext;
		Dataseries* ds,*dsnext;
		Gp* gp,*gpnext;
		gp=firstwindow;
		while (gp) {
			ds=gp->pfirst;
			while (ds) {
				d=ds->pfirst;
				while (d) {
					dnext=d->pnext;
					delete d;
					d=dnext;
				}
				dsnext=ds->pnext;
				delete ds;
				ds=dsnext;
			}
			gpnext=gp->pnext;
			delete gp;
			gp=gpnext;
		}
		nwindow=0;
		firstwindow=NULL;
	}
};


///////////////////////////////////////////////////////////////////////////////////////
// X WINDOWS CLASS DEFINITIONS


// Colour palette for XClient object

const int NCOLOUR=16; // Number of colours in palette

class Colour {

public:
	unsigned long bgr;
	char name[16];
	int id;
	long pixel;

	int blue() {
		return (bgr>>16)/255.0*65535.0;
	}
	int green() {
		return ((bgr&0xFF00)>>8)/255.0*65535.0;
	}
	int red() {
		return (bgr&0xFF)/255.0*65535.0;
	}
};

Colour palette[NCOLOUR] = {
	0x000000,"black",0,0,
	0xFF0000,"dark blue",1,0,
	0x0000FF,"dark red",2,0,
	0x00FF00,"dark green",3,0,
	0xFFFF00,"cyan",4,0,
	0xFF00FF,"violet",5,0,
	0x00FFFF,"yellow",6,0,
	0xAFAFAF,"light grey",7,0,
	0xFF8F8F,"light blue",8,0,
	0x8F8FFF,"light red",9,0,
	0x8FFF8F,"light green",10,0,
	0xFFFF8F,"light cyan",11,0,
	0xFF8FFF,"light violet",12,0,
	0x008F8F,"dark yellow",13,0,
	0x3F3F3F,"dark grey",14,0,
	0xFFFFFF,"white",15,0
};


class XClient; // (forward declaration)

class XPalette {

	friend class XClient;

public:
	Colour* m_palette;
	int m_ncolour;
	Colormap m_colourmap;
	XClient* m_pparent;

	void init(XClient* pparent,Colour* palette,int ncolour);
	long pixel(char* name);
	long pixel(int id);
};


// Menu class for X windows

class XWindow; // (forward declaration)

class XMenu {

public:
	XWindow* m_pparent;
	bool m_visible;
	int x,y,width,height,rowht;
	int nrow;
	char* pstrings;

	XMenu();
	~XMenu();
	void clearmenu();
	void calcsize();
	void addstring(char* text);
	void showmenu(bool ifshow,int ax,int ay);
	void drawmenu();
	int getmenuitem(int ax,int ay);
};


// Individual X window class

// Forward declarations
class Gp;
class CGraph;

class XWindow {

	friend class XClient;

public:
	XClient* m_pparent;
	XWindow* pprev;
	XWindow* pnext;
	unsigned int xpos,ypos;
	unsigned int width,height;
	GC m_gc;
	xtring m_name;
	xtring m_filename;
	Gp* m_pgp; // pointer to corresponding Gp object
	XMenu m_menu;
	bool m_frozen;
	bool m_shown;
	Window m_window;

	XWindow(XClient* pparent,int x,int y,unsigned int w,unsigned int h,char* name);
	void set_window_name();
	void resize_window(unsigned int w,unsigned int h);
	void get_window_specs(int& x,int& y,unsigned int& w,unsigned int& h);
	void copy_bitmap();
	void MoveTo(unsigned int x,unsigned int y);
	void LineTo(unsigned int x,unsigned int y);
	void FillRect(unsigned int x,unsigned int y,unsigned int w,unsigned int h);
	void SetColour(int id);
	void GetTextMetrics(char* text,int& width,int& height);
	int GetFontHeight();
	void TextOut(int xpos,int ypos,char* text);
	void update_window();
	void draw_window(CGraph* pcgraph,bool redraw);
};


// X window manager class

class XClient {

public:// BLARP temporary

	Display* m_display;
	bool m_connected; // true if connection to X server established
	XPalette m_palette;
	int m_ncolour; // number of colours in palette
	XWindow* pfirst; // pointer to start of linked list of XWindow objects
	XWindow* plast;
	Font m_font; // font (NULL for default font)

	void getevent(XEvent& event) {
		XNextEvent(m_display,&event);
	};

	XClient();
	~XClient();
	void disconnect();
	void initialise(char* display_name);
	void save_window_specs(char* filename);
	void retrieve_window_specs(char* filename,char* name,int& x,int& y,
		unsigned int& width,unsigned int& height);
	XWindow* new_window(int x,int y,unsigned int width,unsigned int height,char* name);
};


// Definitions of X class member functions

void XPalette::init(XClient* pparent,Colour* palette,int ncolour) {
	
	// Populate XPalette object with X library-compatible colours

	int c;
	XColor xcolour;
	m_palette=palette;
	m_ncolour=ncolour;
	m_pparent=pparent;
	m_colourmap=DefaultColormap(pparent->m_display,0);
	for (c=0;c<ncolour;c++) {
		xcolour.red=palette[c].red();
		xcolour.green=palette[c].green();
		xcolour.blue=palette[c].blue();
		xcolour.flags=DoRed|DoGreen|DoBlue;
		if (XAllocColor(pparent->m_display,m_colourmap,&xcolour))
			palette[c].pixel=xcolour.pixel;
	}
};

long XPalette::pixel(char* name) {

	// Returns X-library pixel index based on colour name

	int c;
	for (c=0;c<m_ncolour;c++) {
		if (!strcmp(m_palette[c].name,name)) {
			return m_palette[c].pixel;
		}
	}
	return 0; // not found
};

long XPalette::pixel(int id) {

	// Returns X-library pixel index based on palette index
	// 0=black, 15=white

	return m_palette[id].pixel;
}

XMenu::XMenu() {
	pstrings=new char[1];
	if (pstrings) *pstrings='\0';
	m_visible=false;
	nrow=0;
}

XMenu::~XMenu() {

	if (pstrings) delete pstrings;
}

void XMenu::clearmenu() {
	if (pstrings) delete pstrings;
	nrow=0;
	m_visible=false;
	pstrings=NULL;
}

void XMenu::calcsize() {

	int r,w;
	char* pstr=pstrings,*pend;
	width=height=0;
	for (r=0;r<nrow;r++) {
		pend=strstr(pstr,"\n");
		*pend='\0';
		m_pparent->GetTextMetrics(pstr,w,rowht);
		if (w>width) width=w;
		height+=rowht;
		*pend='\n';
		pstr=++pend;
	}
	width+=8;
	height+=4+nrow*4;
}

void XMenu::addstring(char* text) {

// BLARP temporary
//return;

	if (!pstrings) return;
	char* pbuf;
	pbuf=new char[strlen(pstrings)+strlen(text)+2];
	if (pbuf) {
		strcpy(pbuf,pstrings);
		strcat(pbuf,text);
		strcat(pbuf,"\n");
		delete pstrings;
		pstrings=pbuf;
		nrow++;
		calcsize();
	}
}

void XMenu::showmenu(bool ifshow,int ax,int ay) {

	if (ifshow && nrow) {
		if (ax+width<m_pparent->width) x=ax;
		else x=max(m_pparent->width-width,0);
		if (ay+height<m_pparent->height) y=ay;
		else y=max(m_pparent->height-height,0);
		m_visible=true;
	}
	else m_visible=false;
}

void XMenu::drawmenu() {

	if (m_visible) {
		int r,w,h,vpos;
		char* pstr=pstrings,*pend;
		m_pparent->SetColour(7); // light grey background
		m_pparent->FillRect(x,y,width,height);
		m_pparent->SetColour(15);
		m_pparent->MoveTo(x+width,y+1);
		m_pparent->LineTo(x+1,y+1);
		m_pparent->LineTo(x+1,y+height);
		m_pparent->SetColour(14);
		m_pparent->MoveTo(x,y+height);
		m_pparent->LineTo(x+width,y+height);
		m_pparent->LineTo(x+width,y);
		m_pparent->SetColour(0); // black text
		vpos=y;
		for (r=0;r<nrow;r++) {
			pend=strstr(pstr,"\n");
			*pend='\0';
			vpos+=rowht+4;
			m_pparent->TextOut(x+5,vpos,pstr);
			*pend='\n';
			pstr=++pend;
		}		
	}
}

int XMenu::getmenuitem(int ax,int ay) {

	if (m_visible) {
		if (ax>=x && ax<x+width && ay>=y+2 && ay<y+height-2)
			return (ay-(y+2))/(rowht+4);
		else
			return -1; // not found
	}
	else return -1;
}


XWindow::XWindow(XClient* pparent,int x,int y,unsigned int w,unsigned int h,char* name) {

	// Constructor for new X window

	XSetWindowAttributes att;
	XGCValues gcatt;
	m_name=name;
	xtring ncopy=name;
	char* pchar=(char*)ncopy;
	while (*pchar) {
		if (*pchar==' ') strcpy(pchar,pchar+1);
		else pchar++;
	}
	m_filename=ncopy;
	width=w;
	height=h;
	m_frozen=false;
	att.event_mask=ExposureMask|StructureNotifyMask|ButtonPressMask;
	att.colormap=pparent->m_palette.m_colourmap;
	att.background_pixel=pparent->m_palette.pixel("white");
	att.cursor=XCreateFontCursor(pparent->m_display,68);
	m_pparent=pparent;
	m_window=XCreateWindow(pparent->m_display,
		DefaultRootWindow(pparent->m_display),
		x,y,width,height,0,
		CopyFromParent,InputOutput,CopyFromParent,
		CWEventMask|CWBackPixel|CWColormap|CWCursor,&att);
	pprev=pparent->plast;
	pnext=NULL;
	pparent->plast=this;
	if (!pparent->pfirst) pparent->pfirst=this;
	else pprev->pnext=this;
	gcatt.background=pparent->m_palette.pixel("white");
	gcatt.foreground=pparent->m_palette.pixel("black");
	m_gc=XCreateGC(pparent->m_display,m_window,
		GCBackground|GCForeground,&gcatt);
	XFlushGC(pparent->m_display,m_gc);
	if (pparent->m_font) XSetFont(pparent->m_display,m_gc,pparent->m_font);
	m_menu.m_pparent=this;
	m_shown=false;
}

void XWindow::set_window_name() {
	const unsigned char* pname=(const unsigned char*)(char*)m_name;
	XChangeProperty(m_pparent->m_display,m_window,XA_WM_NAME,XA_STRING,8,
		PropModeReplace,pname,m_name.len());	
}


void XWindow::resize_window(unsigned int w,unsigned int h) {
	width=w;
	height=h;
}

void XWindow::get_window_specs(int& x,int& y,unsigned int& w,unsigned int& h) {

	XWindowAttributes att;
	XGetWindowAttributes(m_pparent->m_display,m_window,&att);
	if (att.map_state==IsViewable) {
		x=att.x;
		y=att.y;
		w=att.width;
		h=att.height;
	}
	else {
		x=0;
		y=0;
		w=width;
		h=height;
	}
}

void XWindow::copy_bitmap() {

	XWindowAttributes att;
	XGetWindowAttributes(m_pparent->m_display,m_window,&att);
	Pixmap p=XCreatePixmap(m_pparent->m_display,m_window,width,height,att.depth);
	XWriteBitmapFile(m_pparent->m_display,"window.bmp",p,width,height,-1,-1);
	XFreePixmap(m_pparent->m_display,p);
}

void XWindow::update_window() {

	// Called to inform XClient object that window needs updating / drawing

	if (m_frozen) return;
	XEvent event;
	XExposeEvent expose_event;
	expose_event.type=Expose;
	expose_event.window=m_window;
	event.xexpose=expose_event;
	XSendEvent(m_pparent->m_display,m_window,False,ExposureMask,&event);
	XFlush(m_pparent->m_display);
};

void XWindow::MoveTo(unsigned int x,unsigned int y) {

	// Starting position in pixels for next line draw
	// (0,0) is top-left corner

	xpos=x;
	ypos=y;
};

void XWindow::LineTo(unsigned int x,unsigned int y) {

	// Line from current position (set by MoveTo) to specified point

	XDrawLine(m_pparent->m_display,m_window,m_gc,xpos,ypos,x,y);
	MoveTo(x,y);
};

void XWindow::FillRect(unsigned int x,unsigned int y,unsigned int w,unsigned int h) {

	// Fill rectangle with foreground colour

	XFillRectangle(m_pparent->m_display,m_window,m_gc,x,y,w,h);
};

void XWindow::SetColour(int id) {

	// Sets foreground colour based on palette colour index

	XSetForeground(m_pparent->m_display,m_gc,m_pparent->m_palette.pixel(id));
};

int XWindow::GetFontHeight() {

	// Returns window font height in pixels from baseline

	XFontStruct* pfstruc;
	XGCValues gcv;
	int ascent,descent;
	XCharStruct box;
	XGetGCValues(m_pparent->m_display,m_gc,GCFont,&gcv);
	pfstruc=XQueryFont(m_pparent->m_display,gcv.font);
	if (pfstruc) {
		ascent=pfstruc->ascent;
		descent=pfstruc->descent;
	}
	else ascent=descent=NULL;
	return ascent;
};

void XWindow::GetTextMetrics(char* text,int& width,int& height) {

	// Returns width and height of specified text in pixels

	int dir,ascent,descent;
	XCharStruct box;
	XQueryTextExtents(m_pparent->m_display,XGContextFromGC(m_gc),text,strlen(text),&dir,
		&ascent,&descent,&box);
	width=box.width;
	height=ascent*1.4;
};

void XWindow::TextOut(int xpos,int ypos,char* text) {

	// Print text from specified pixel position
	// Anchor point is baseline just left of first character in text

	int dir,ascent,descent;
	XCharStruct box;
	XQueryTextExtents(m_pparent->m_display,XGContextFromGC(m_gc),text,strlen(text),&dir,
		&ascent,&descent,&box);
	
	XTextItem t;
	t.chars=text;
	t.nchars=strlen(text);
	t.delta=0;
	t.font=None;
	XDrawText(m_pparent->m_display,m_window,m_gc,xpos,ypos-ascent*0.2,&t,1);
};

void XWindow::draw_window(CGraph* pcgraph,bool redraw) {

	// Code for drawing the contents of this window
	// Should be called by XClient object

	if (m_pgp->m_update || redraw || ifmono) {
		m_pgp->m_update=true;
		XClearWindow(m_pparent->m_display,m_window);
	}
	pcgraph->update(m_pgp);

	m_menu.drawmenu();

	XMapWindow(m_pparent->m_display,m_window);
	XFlush(m_pparent->m_display);
};

XClient::XClient() {

	// Constructor for X windows manager object

	m_connected=false;
	pfirst=plast=NULL;
};

XClient::~XClient() {

	// Deconstructor for X windows manager object

	disconnect();
	XWindow* pthis,*pnext;
	pthis=pfirst;
	while (pthis) {
		pnext=pthis->pnext;
		delete pthis;
		pthis=pnext;
	}
};

void XClient::disconnect() {
	
	// Will switch off thread function listen

	m_connected=false;
}


void XClient::initialise(char* display_name) {

	// Initialise X windows manager
	// Establishes connection to specified X windows server

	XFontStruct *pfstruc;
	m_display=XOpenDisplay(display_name);
	if (m_display) {
		m_connected=true;
		m_ncolour=NCOLOUR;
		m_palette.init((XClient*)this,palette,m_ncolour);

		// Load font and get metrics (try to get Courier or similar)

		pfstruc=XLoadQueryFont(m_display,"*courier-medium-r-normal*-10-*");
		if (pfstruc) m_font=pfstruc->fid;
		else {
			m_font=NULL;
		}
	}
};

XWindow* XClient::new_window(int x,int y,unsigned int width,unsigned int height,char* name) {

	// Should be called by application to create a new X window

	if (m_connected) {
		return new XWindow((XClient*)this,x,y,width,height,name);
	}
	// else
	return NULL;
};

void XClient::save_window_specs(char* filename) {

	// Saves name,x,y,width,height for all windows to a binary file

	XWindow* pthis;
	unsigned int nwindow=0;
	unsigned int length;
	int x,y;
	unsigned int w,h;
	char flag='0';

	if (m_connected && pfirst) {
		
		FILE* out=fopen(filename,"wb");
		if (out) {
			pthis=pfirst;
			while (pthis) {
				nwindow++;
				pthis=pthis->pnext;
			}
			fwrite(&flag,sizeof(char),1,out);
			fwrite(&nwindow,sizeof(unsigned int),1,out);
			pthis=pfirst;
			while (pthis) {
				length=pthis->m_name.len();
				fwrite(&length,sizeof(unsigned int),1,out);
				fwrite((char*)(pthis->m_name),sizeof(char),length+1,out);
				pthis->get_window_specs(x,y,w,h);
				fwrite(&x,sizeof(int),1,out);
				fwrite(&y,sizeof(int),1,out);
				fwrite(&w,sizeof(unsigned int),1,out);
				fwrite(&h,sizeof(unsigned int),1,out);
				pthis=pthis->pnext;
			}
			rewind(out);
			flag='1';
			fwrite(&flag,sizeof(char),1,out);
			fclose(out);
		}
	}
}

void XClient::retrieve_window_specs(char* filename,char* name,int& x,int& y,
	unsigned int& width,unsigned int& height) {

	unsigned int n,siz,w,h;
	int i,ax,ay;
	xtring buf;
	char flag;

	FILE* in=fopen(filename,"rb");
	if (in) {
		fread(&flag,sizeof(char),1,in);
		if (flag!='1') {
			fclose(in);
			return;
		}
		fread(&n,sizeof(unsigned int),1,in);
		for (i=0;i<n;i++) {
			fread(&siz,sizeof(unsigned int),1,in);
			buf.reserve(siz);
			fread((char*)buf,sizeof(char),siz+1,in);
			fread(&ax,sizeof(int),1,in);
			fread(&ay,sizeof(int),1,in);
			fread(&w,sizeof(unsigned int),1,in);
			fread(&h,sizeof(unsigned int),1,in);
			if (buf==name) {
				x=ax;
				y=ay;
				width=w;
				height=h;
				fclose(in);
				return;
			}
		}
		fclose(in);
	}
}


///////////////////////////////////////////////////////////////////////////////////////
// X WINDOWS COMPATIBLE GRAPHICS LIBRARY DEFINITIONS


void plotlinesection(XWindow* pwin,int symbol,int fromx,int fromy,int tox,int toy,
	double& progress,double scale) {

	// Plots a section of a monochrome line pattern
	// (Global function)

	int i,section;
	double angle,length,spointer,sum,slength,rlength,tlength;
	double dfromx=fromx,dfromy=fromy,dtox=tox,dtoy=toy,ddtox,ddtoy;
	double rscaled;
	bool found,slopy;

	scale*=SCALE_LINE_SYMBOL;

	length=0.0;
	for (i=0;i<linesymbol[symbol].nsection;i++) length+=linesymbol[symbol].section[i];

	if (fromx==tox)
		slopy=false;
	else
		slopy=true;

	if (slopy) {
		angle=atan((dtoy-dfromy)/(dtox-dfromx));
		tlength=sqrt((dtoy-dfromy)*(dtoy-dfromy)+(dtox-dfromx)*(dtox-dfromx));
	}
	else
		tlength=dtoy-dfromy;
	section=0;
	spointer=0.0;
	sum=0.0;
	found=false;
	while (!found) {
		slength=(double)linesymbol[symbol].section[section]/length;
		if (sum+slength>=progress) {
			found=true;
			rlength=sum+slength-progress;
		}
		else {
			section++;
			sum+=slength;
		}
	}
	while (true) {
		rscaled=rlength*length*scale;
		if (rscaled<tlength) { // finish plotting this section and move on to next section
			if (slopy) {
				ddtox=dfromx+rscaled*cos(angle);
				ddtoy=dfromy+rscaled*sin(angle);
			}
			else {
				ddtox=dfromx;
				ddtoy=dfromy+rscaled;
			}
			if (!(section%2)) { // even section - pen down
				pwin->MoveTo(dfromx+0.5,dfromy+0.5);
				pwin->LineTo(ddtox+0.5,ddtoy+0.5);
			}
			progress+=rlength;
			tlength-=rscaled;
			section++;
			if (section>linesymbol[symbol].nsection) {
				section=0;
				progress=0.0;
				sum=0.0;
			}
			dfromx=ddtox;
			dfromy=ddtoy;
			rlength=(double)linesymbol[symbol].section[section]/length;
		}
		else { // plot this section to end and return
			if (rlength) {
				if (!(section%2)) { // even section - pen down
					pwin->MoveTo(dfromx+0.5,dfromy+0.5);
					pwin->LineTo(dtox+0.5,dtoy+0.5);
				}
				progress+=rlength*(tlength/rlength/length/scale);
			}
			return;
		}
	}
}

double sign(double d) {

	// Returns sign (+/- 1) of d

	if (!d) return 1.0;
	return d/fabs(d);
}

void CGraph::tag(int n) {

	// Adds a tag to all currently untagged windows

	xtring name;

	if (nwindow) {
		Gp* pgp=firstwindow;
		if (!pgp->tag) {
			pgp->tag=n;
			name.printf("(%s from run %d)",pgp->name,n);
			strcpy(pgp->name,(char*)name);
		}
		while (pgp->pnext) {
			pgp=pgp->pnext;
			if (!pgp->tag) {
				pgp->tag=n;
				name.printf("(%s from run %d)",pgp->name,n);
				strcpy(pgp->name,(char*)name);
			}
		}
	}
}

bool CGraph::createwindow(XClient* pxclient,char* wname,unsigned int width,
	unsigned int height,bool use_saved_dimensions) {

	// Creates a new output window with specified name, width and height

	int i;

	if (nwindow) { // return if graph with this name already exists
		xtring n1,n2;
		n1=wname;
		n1=n1.upper();
		Gp* pgp=firstwindow;
		n2=pgp->name;
		n2=n2.upper();
		if (n1==n2) return true;
		while (pgp->pnext) {
			pgp=pgp->pnext;
			n1=wname;
			n1=n1.upper();
			n2=pgp->name;
			n2=n2.upper();
			if (n1==n2) return true;
		}
	}

	Gp* newgp=new Gp(wname,width,height,0);
	if (!newgp) {
		return false;
	}
	newgp->pnext=NULL;

	if (!nwindow) {
		firstwindow=newgp;
	}
	else {
		Gp* pgp=firstwindow;
		for (i=0;i<nwindow;i++) {
			if (!pgp->pnext) pgp->pnext=newgp;
			else pgp=pgp->pnext;
		}
	}

	nwindow++;

	// Set parameters to send to XClient object
	
	int x=0,y=0;
	if (use_saved_dimensions)
		pxclient->retrieve_window_specs(WSFILE,wname,x,y,width,height);
	XWindow* pnewwnd=pxclient->new_window(x,y,width,height,wname);
	if (pnewwnd) {
		newgp->pwnd=pnewwnd;
		pnewwnd->m_pgp=newgp;
	}
	else {
		delete newgp;
		return false;
	}

	pnewwnd->update_window();
	
	newgp->m_update=true;
	
	return true;
}

Gp& CGraph::getgp(char* name) {

	// Returns reference to Gp object for graph with given name

	xtring n1,n2;

	n1=name;
	n1=n1.upper();

	Gp* pgp=firstwindow;
	while (pgp) {
		n2=pgp->name;
		n2=n2.upper();
		if (n1==n2) return *pgp;
		pgp=pgp->pnext;
	}
	return *pgp;
}

void CGraph::setlabelstyle(Labeldef& labelstyle,int width,int places,bool exp) {

	// Sets x-axis label style parameters

	labelstyle.width=width;
	labelstyle.places=places;
	labelstyle.exp=exp;
	return;

}

void CGraph::autoscaleaxis(Axisdef& axispars,Labeldef& labelstyle,double& base,double& mul,
		double& setmin,double& setmax,double& autotic,double& minortic,
		double min,double max,double round,bool fixrange) {

	// Automatically scales graph area in x direction
	// if zerox non-zero (true) and minx>=0, x origin set to 0
	// if round<10 used as a criterion for rounding (e.g., if
	// round=2, rounding is to nearest 0.002, 0.2, 2, 20 etc.
	// if round>=10, an internal rule is used instead

	if (axispars.zero && !axispars.log) {
		if (min<=0.0 && max<=0.0) max=0.0;
		else if (min>=0.0 && max>=0.0) min=0.0;
	}

	double unit,gmax,gmin,minorunit;
	double interval=max-min;
	double mord;

	if (!interval) {
		if (!max) {
			interval=0.001;
			max=min+0.001;
		}
		else {
			interval=0.001;
			min=max-0.001;
		}
	}

	if (axispars.log) {
		autotic=1;
		if (interval>=MAXBGAUTOTIC) minortic=0.0;
		else minortic=1.0; // between log 1.0 and log 10.0
	}
	else {
		mord=floor(log10(interval));
		interval/=pow(10.0,mord);	// now in the range 1.0-9.9999

		if (round>=10.0) {
			unit=0.2;
			if (interval/unit>MAXBGAUTOTIC) unit=0.5;
			if (interval/unit>MAXBGAUTOTIC) unit=1.0;
			if (interval/unit>MAXBGAUTOTIC) unit=2.0;
		}
		else {
			unit=round;
			if (interval/unit<5.0) unit/=10.0;
		}

		if (unit==0.3 || unit==0.5) minorunit=0.1;
		else if (unit==3.0 || unit==5.0) minorunit=1.0;
		else if (unit==1.0) minorunit=0.2;
		else minorunit=unit/2.0;

		// now unit*10^mord is the automatic tic interval

		autotic=unit*pow(10.0,mord);
		minortic=minorunit*pow(10.0,mord);
	}

	if (fixrange) {
		gmax=max;
		gmin=min;
	}
	else {
		gmax=floor(max/autotic+0.999999)*autotic;
		gmin=floor(min/autotic+0.000001)*autotic;
	}

	// Set axis parameters
	base=gmin;
	mul=1.0/(gmax-gmin);
	setmin=gmin;
	setmax=gmax;

	if (autotic>999999.0) {
		setlabelstyle(labelstyle,0,1,true);
		return;
	}
	if (autotic>=1.0) {
		setlabelstyle(labelstyle,0,0,false);
		return;
	}
	if (autotic>=0.09999) {
		setlabelstyle(labelstyle,0,1,false);
		return;
	}
	if (autotic>=0.009999) {
		setlabelstyle(labelstyle,0,2,false);
		return;
	}
	if (autotic>=0.0009999) {
		setlabelstyle(labelstyle,0,3,false);
		return;
	}
	if (autotic>=0.00009999) {
		setlabelstyle(labelstyle,0,4,false);
		return;
	}
	setlabelstyle(labelstyle,0,1,true);
}


void CGraph::scalex(char* name,double min,double max,double round) {

	// Scales x-axis on specified graph window from approximately min to max
	// with specified rounding

	Gp& gp=getgp(name);
	
	double xbase_old=gp.xbase;
	double xmul_old=gp.xmul;
	double xmin_old=gp.xmin;
	double xmax_old=gp.xmax;
	
	autoscaleaxis(gp.xaxispars,gp.xlabelstyle,gp.xbase,gp.xmul,gp.xmin,gp.xmax,
		gp.autoticx,gp.minorticx,min,max,round,false);
		
	if (gp.xbase!=xbase_old || gp.xmul!=xmul_old || gp.xmin!=xmin_old || gp.xmax!=xmax_old)
		gp.m_update=true;
}

void CGraph::scaley(char* name,double min,double max,double round) {

	// Scales y-axis on specified graph window from approximately min to max
	// with specified rounding

	Gp& gp=getgp(name);
	
	double ybase_old=gp.ybase;
	double ymul_old=gp.ymul;
	double ymin_old=gp.ymin;
	double ymax_old=gp.ymax;
	
	autoscaleaxis(gp.yaxispars,gp.ylabelstyle,gp.ybase,gp.ymul,gp.ymin,gp.ymax,
		gp.autoticy,gp.minorticy,min,max,round,false);

	if (gp.ybase!=ybase_old || gp.ymul!=ymul_old || gp.ymin!=ymin_old || gp.ymax!=ymax_old)
		gp.m_update=true;
}

void CGraph::logprint(char* labtext,double ctic,bool exp) {

	// Special print formats for log-axis labels

	double label=pow(10.0,ctic);
	if (ctic>5.9999 || ctic<=-4.0001 || exp) bgsprint(labtext,label,0,0,1);
	else if (ctic>-0.00001) bgsprint(labtext,label,0,0,0);
	else if (ctic>-1.0001) bgsprint(labtext,label,0,1,0);
	else if (ctic>-2.0001) bgsprint(labtext,label,0,2,0);
	else if (ctic>-3.0001) bgsprint(labtext,label,0,3,0);
	else if (ctic>-4.0001) bgsprint(labtext,label,0,4,0);
}

void CGraph::bgsprint(char* string,double var,int width,int places,bool exp) {

	// Prints "var" to "string" with specified width, places and exp

	char buf[50],buf2[15],*bufptr,*strptr;
	int flag;

	if (exp) sprintf(buf2,"%%%d.%dE",width,places);
	else sprintf(buf2,"%%%d.%df",width,places);
	sprintf(buf,buf2,var);

	// check if "-0.000E+12" or something silly like that

	bufptr=buf;
	while (*bufptr==' ') bufptr++;
	if (*bufptr=='-') { // check if only zero digits to the E
		flag=0;
		while (*bufptr!='E' && *bufptr) {
			if (*bufptr>='1' && *bufptr<='9') flag=1;  //a non-zero digit
			bufptr++;
		}
		if (!flag) {
			bgsprint(string,0.0,width,places,exp);  //rerun with var=exactly 0
			return;
		}
	}
	strcpy(string,buf);
}


void CGraph::pline(XWindow* pwin,int hstart,int vstart,int hend,int vend) {

	// Draws a line from (hstart,vstart) to (hend,vend) on specified X window

	pwin->MoveTo(hstart,vstart);
	pwin->LineTo(hend,vend);
	pwin->LineTo(hstart,vstart);
}

void CGraph::ptext(XWindow* pwin,Gp& gp,int hpos,int vpos,justtype just,char* text) {

	// Prints specifed text justified as specified relative to anchor point (hpos,vpos)
	// on specified X window

	int cx,cy;
	pwin->GetTextMetrics(text,cx,cy);

	if (just==AUTOJUST || just==CENTRE || just==CC)
		pwin->TextOut(hpos-cx/2,vpos+cy/2*gp.dir,text);
	else if (just==LEFT || just==CL)
		pwin->TextOut(hpos,vpos+cy/2*gp.dir,text);
	else if (just==RIGHT || just==CR)
		pwin->TextOut(hpos-cx,vpos+cy/2*gp.dir,text);
	else if (just==TL)
		pwin->TextOut(hpos,vpos+cy*gp.dir,text);
	else if (just==TC)
		pwin->TextOut(hpos-cx/2,vpos+cy*gp.dir,text);
	else if (just==TR)
		pwin->TextOut(hpos-cx,vpos+cy*gp.dir,text);
	else if (just==BL)
		pwin->TextOut(hpos,vpos,text);
	else if (just==BC)
		pwin->TextOut(hpos-cx/2,vpos,text);
	else // just==BR
		pwin->TextOut(hpos-cx,vpos,text);
}

void CGraph::plotdata(XWindow* pwin,Gp& gp,Dataseries& series) {

	// Outputs a curve on a specified X window

	Data* pdata=series.pfirst;

	if (!pdata) return;

	if (!ifmono) {

		pwin->SetColour(series.colourid);

		pwin->MoveTo(gp.gtopx(pdata->x),gp.gtopy(pdata->y));
		while (pdata->pnext) {
			pdata=pdata->pnext;
			if (pdata->plotted && !gp.m_update) {
				if (pdata->pnext) {
					if (!pdata->pnext->plotted)
						pwin->MoveTo(gp.gtopx(pdata->x),gp.gtopy(pdata->y));
				}
			}
			else
				pwin->LineTo(gp.gtopx(pdata->x),gp.gtopy(pdata->y));
			pdata->plotted=true;
		}
	}
	else {

		double progress=0.0;
		int fromx,fromy,tox,toy;

		pwin->SetColour(0); // black

		fromx=gp.gtopx(pdata->x);
		fromy=gp.gtopy(pdata->y);
		while (pdata->pnext) {
			pdata=pdata->pnext;
			tox=gp.gtopx(pdata->x);
			toy=gp.gtopy(pdata->y);
			plotlinesection(pwin,series.series%NLINESYMBOL,fromx,fromy,tox,toy,progress,gp.scale);
				// use tic length as scaling factor
			pdata->plotted=true;
			fromx=tox;
			fromy=toy;
		}
	}
}

void CGraph::legend(XWindow* pwin,Gp& gp) {

	// Outputs a legend on a specified X window

	double progress;
	int cx,cy;

	int ht=0;
	int maxwidth=0;
	int textht;
	int avail,needed;
	int tlx,tly,y;

	Dataseries* pseries=gp.pfirst;

	if (!pseries) return;

	while (pseries) {
		pwin->GetTextMetrics(pseries->name,cx,cy);
		if (cx>maxwidth) maxwidth=cx;
		ht+=cy;
		pseries=pseries->pnext;
	}

	textht=cy;

	maxwidth+=gp.leglinelen+textht/2;

	// adjust graph shape if necessary to fit legend;

	if ((double)maxwidth/(double)ht>(double)gp.width/(double)gp.height) { // place under graph

		avail=(gp.margin+gp.height)*gp.dir-gp.vmax-textht-gp.ticstyle.length
			-gp.ticstyle.length/4;
		needed=ht+textht;
		if (avail<needed && gp.vrange>needed-avail+100) {
			gp.vmax-=(needed-avail)*gp.dir;
			gp.vrange-=needed-avail;
		}
		tly=gp.vmax+(textht+textht/2+gp.ticstyle.length+gp.ticstyle.length/4)*gp.dir;
		tlx=gp.margin+gp.width/2-maxwidth/2;
	}
	else { // place beside graph

		avail=gp.width+gp.margin-gp.hmin-gp.hrange;
		needed=maxwidth+2*textht;
		if (avail<needed && gp.hrange>needed-avail+100) gp.hrange-=needed-avail;
		tlx=gp.hmin+gp.hrange+textht;
		tly=(gp.margin+gp.height/2-ht/2)*gp.dir;
	}

	// Now plot legend

	pseries=gp.pfirst;

	y=tly;
	while (pseries) {

		if (ifmono) {

			pwin->SetColour(0); // black
			progress=0.0;
			plotlinesection(pwin,pseries->series%NLINESYMBOL,tlx,y+textht/2*gp.dir,
				tlx+gp.leglinelen,y+textht/2*gp.dir,progress,gp.scale);
		}
		else {
			pwin->SetColour(pseries->colourid);
			pwin->MoveTo(tlx,y+textht/2*gp.dir);
			pwin->LineTo(tlx+gp.leglinelen,y+textht/2*gp.dir);
		}
	
		ptext(pwin,gp,tlx+gp.leglinelen+textht,y,TL,pseries->name);
		y+=textht*gp.dir;
		pseries=pseries->pnext;
	}
}

void CGraph::xaxis(XWindow* pwin,Gp& gp) {
	
	// Outputs an X axis on a specified X window

	pwin->SetColour(0);
	pline(pwin,gp.hmin,gp.vmax,gp.hmin+gp.hrange,gp.vmax);

}

void CGraph::yaxis(XWindow* pwin,Gp& gp) {

	// Outputs a Y axis on a specified X window

	pwin->SetColour(0);
	pline(pwin,gp.hmin,gp.vmax,gp.hmin,gp.vmax-gp.vrange);
}

void CGraph::xtics(XWindow* pwin,Gp& gp) {

	// Plots tics and labels, with xdisplacement (graph units) if desired
	// Displaced tics are shown only if they fall within the axis length

	double gmin,gmax,tol,ctic,axend,hh,cmtic;
	Location s,f,sm,fm;
	double yintercept;
	int ydispl;
	int mult=1;
	char labtext[50];

	pwin->SetColour(0);

	yintercept=gp.ptogy(gp.vmax);
	ydispl=gp.xaxispars.displ;

	gmin=gp.ptogx(gp.hmin);
	gmax=gp.ptogx(gp.hmin+gp.hrange);
	tol=(gmax-gmin)*0.001;
	if ((double)(gp.vmax-ydispl-gp.gtopy(yintercept))>(double)gp.vrange/2.0) mult=-1;
	if (gp.ticstyle.centre) {
		s.v=gp.gtopy(yintercept)+gp.ticstyle.length/2+ydispl;
		f.v=gp.gtopy(yintercept)-gp.ticstyle.length/2+ydispl;
		sm.v=gp.gtopy(yintercept)+gp.minorticstyle.length/2+ydispl;
		fm.v=gp.gtopy(yintercept)-gp.minorticstyle.length/2+ydispl;
	}
	else {
		s.v=gp.gtopy(yintercept)+ydispl;
		f.v=gp.gtopy(yintercept)+ydispl+gp.ticstyle.length*mult;
		sm.v=s.v;
		fm.v=gp.gtopy(yintercept)+ydispl+gp.minorticstyle.length*mult;
	}

	axend=gmin+tol*sign(gmin);
	ctic=(double)long(axend/gp.autoticx)*gp.autoticx;
	while (ctic<gmax+tol) {
		if (ctic>gmin-tol) {
			s.h=f.h=gp.gtopx(ctic);
			pline(pwin,s.h,s.v,f.h,f.v);
			
			if (gp.xaxispars.log) logprint(labtext,ctic,gp.xlabelstyle.exp);
			else bgsprint(labtext,ctic,gp.xlabelstyle.width,
				gp.xlabelstyle.places,gp.xlabelstyle.exp);

			ptext(pwin,gp,f.h,f.v+gp.ticstyle.length/4,TC,labtext);
			if (gp.wantminortics) {
				if (gp.xaxispars.log) {
					hh=gp.minorticx;
					cmtic=ctic+log10(hh);
					while (cmtic<ctic+gp.autoticx-tol && cmtic<gmax) {
						sm.h=fm.h=gp.gtopx(cmtic);
						pline(pwin,s.h,s.v,f.h,f.v);
						hh+=gp.minorticx;
						cmtic=ctic+log10(hh);
					}
				}
				else {
					cmtic=ctic+gp.minorticx;
					while (cmtic<ctic+gp.autoticx-tol && cmtic<gmax) {
						sm.h=fm.h=gp.gtopx(cmtic);
						pline(pwin,sm.h,sm.v,fm.h,fm.v);
						cmtic+=gp.minorticx;
					}
				}
			}
		}
		ctic+=gp.autoticx;
	}
}

void CGraph::ytics(XWindow* pwin,Gp& gp) {

	// Plots tics and labels, with xdisplacement (graph units) if desired
	// Displaced tics are shown only if they fall within the axis length

	double gmin,gmax,tol,ctic,axend,hh,cmtic;
	Location s,f,sm,fm;
	double xintercept;
	int xdispl;
	int mult=1;
	char labtext[50];

	pwin->SetColour(0);

	xintercept=gp.ptogx(gp.hmin);
	xdispl=gp.yaxispars.displ;

	gmin=gp.ptogy(gp.vmax);
	gmax=gp.ptogy(gp.vmax-gp.vrange);
	tol=(gmax-gmin)*0.001;
	if ((double)(gp.gtopx(xintercept)-gp.hmin+xdispl)>(double)gp.hrange/2.0) mult=-1;
	if (gp.ticstyle.centre) {
		s.h=gp.gtopx(xintercept)+gp.ticstyle.length/2+xdispl;
		f.h=gp.gtopx(xintercept)-gp.ticstyle.length/2+xdispl;
		sm.h=gp.gtopx(xintercept)+gp.minorticstyle.length/2+xdispl;
		fm.h=gp.gtopx(xintercept)-gp.minorticstyle.length/2+xdispl;
	}
	else {
		s.h=gp.gtopx(xintercept)-xdispl;
		f.h=gp.gtopx(xintercept)-xdispl-gp.ticstyle.length*mult;
		sm.h=s.h;
		fm.h=gp.gtopx(xintercept)-xdispl-gp.minorticstyle.length*mult;
	}

	axend=gmin+tol*sign(gmin);
	ctic=(double)long(axend/gp.autoticy)*gp.autoticy;
	while (ctic<gmax+tol) {
		if (ctic>gmin-tol) {
			s.v=f.v=gp.gtopy(ctic);
			pline(pwin,s.h,s.v,f.h,f.v);
			
			if (gp.yaxispars.log) logprint(labtext,ctic,gp.ylabelstyle.exp);
			else bgsprint(labtext,ctic,gp.ylabelstyle.width,
				gp.ylabelstyle.places,gp.ylabelstyle.exp);

			ptext(pwin,gp,f.h-gp.ticstyle.length/4,f.v,CR,labtext);
			if (gp.wantminortics) {
				if (gp.yaxispars.log) {
					hh=gp.minorticy;
					cmtic=ctic+log10(hh);
					while (cmtic<ctic+gp.autoticy-tol && cmtic<gmax) {
						sm.v=fm.v=gp.gtopy(cmtic);
						pline(pwin,s.h,s.v,f.h,f.v);
						hh+=gp.minorticy;
						cmtic=ctic+log10(hh);
					}
				}
				else {
					cmtic=ctic+gp.minorticy;
					while (cmtic<ctic+gp.autoticy-tol && cmtic<gmax) {
						sm.v=fm.v=gp.gtopy(cmtic);
						pline(pwin,sm.h,sm.v,fm.h,fm.v);
						cmtic+=gp.minorticy;
					}
				}
			}
		}
		ctic+=gp.autoticy;
	}
}

bool CGraph::newseries(char* wname,char* name) {

	// Attempts to start a new series (with specified name) on
	// Graph window with specified wname

	xtring n1,n2;
	int ncolour;

	Gp& gp=getgp(wname);
	ncolour=gp.pwnd->m_pparent->m_ncolour;

	n1=name;
	n1=n1.upper();

	// Check if this series already exists
	Dataseries* pds=gp.pfirst;

	while (pds) {
		n2=pds->name;
		n2=n2.upper();
		if (n1==n2) return true; // return if series already exists
		pds=pds->pnext;
	}

	Dataseries* pnewseries=new Dataseries(gp.nseries,gp.nseries%(ncolour-1),name);
	if (!pnewseries) return false;

	gp.nseries++;

	if (gp.pfirst) {
		pds=gp.pfirst;
		while (pds->pnext) pds=pds->pnext;
		pds->pnext=pnewseries;
	}
	else gp.pfirst=pnewseries;

	gp.m_update=true;
	
	return true;
}

bool CGraph::data(char* wname,char* name,double xval,double yval,bool rescale) {

	// Attempts to add a data coordinate to specified graph window (wname) and
	// data series (name)

	Gp& gp=getgp(wname);
	xtring n1,n2;

	n1=name;
	n1=n1.upper();

	// Try and find series
	Dataseries* pds=gp.pfirst;
	Dataseries* thisseries=NULL;

	while (pds && !thisseries) {
		n2=pds->name;
		n2=n2.upper();
		if (n1==n2) thisseries=pds;
		pds=pds->pnext;
	}

	if (!thisseries) return false;

	Data* newdatum=new Data;
	if (!newdatum) return false;

	newdatum->pnext=NULL;
	newdatum->x=xval;
	newdatum->y=yval;
	newdatum->plotted=false;

	double xmin=xval;
	double xmax=xval;
	double ymin=yval;
	double ymax=yval;

	Data* pdata;
	if (thisseries->pfirst) {
		pdata=thisseries->pfirst;
		if (pdata->x>xmax) xmax=pdata->x;
		if (pdata->x<xmin) xmin=pdata->x;
		if (pdata->y>ymax) ymax=pdata->y;
		if (pdata->y<ymin) ymin=pdata->y;
		while (pdata->pnext) {
			pdata=pdata->pnext;
			if (pdata->x>xmax) xmax=pdata->x;
			if (pdata->x<xmin) xmin=pdata->x;
			if (pdata->y>ymax) ymax=pdata->y;
			if (pdata->y<ymin) ymin=pdata->y;
		}
		pdata->pnext=newdatum;
	}
	else thisseries->pfirst=newdatum;

	if (gp.ndata<2) {
		scalex(wname,xmin,xmax,100.0);
		scaley(wname,ymin,ymax,100.0);
		
		//autoscaleaxis(gp.xaxispars,gp.xlabelstyle,gp.xbase,gp.xmul,gp.xmin,gp.xmax,
		//	gp.autoticx,gp.minorticx,xmin,xmax,100.0,false);
		//autoscaleaxis(gp.yaxispars,gp.ylabelstyle,gp.ybase,gp.ymul,gp.ymin,gp.ymax,
		//	gp.autoticy,gp.minorticy,ymin,ymax,100.0,false);
	}
	else if (rescale) {
		if (xval<gp.xmin || xval>gp.xmax)
			scalex(wname,min(xval,gp.xmin),max(xval,gp.xmax),100.0);
		if (yval<gp.ymin || yval>gp.ymax)
			scaley(wname,min(yval,gp.ymin),max(yval,gp.ymax),100.0);
	
		//if (xval<gp.xmin || xval>gp.xmax)
		//	autoscaleaxis(gp.xaxispars,gp.xlabelstyle,gp.xbase,gp.xmul,gp.xmin,gp.xmax,
		//		gp.autoticx,gp.minorticx,min(xval,gp.xmin),max(xval,gp.xmax),100.0,false);
		//if (yval<gp.ymin || yval>gp.ymax)
		//	autoscaleaxis(gp.yaxispars,gp.ylabelstyle,gp.ybase,gp.ymul,gp.ymin,gp.ymax,
		//		gp.autoticy,gp.minorticy,min(yval,gp.ymin),max(yval,gp.ymax),100.0,false);
	}

	gp.ndata++;

	// call window to update

	gp.pwnd->update_window();

	return true;
}

void CGraph::update(Gp* pgp) {

	// Outputs graph on specified X Window

	Dataseries* pseries;

	Gp& gp=*pgp;

	if (gp.m_update) {
		legend(pgp->pwnd,gp);
		xaxis(pgp->pwnd,gp);
		yaxis(pgp->pwnd,gp);
		xtics(pgp->pwnd,gp);
		ytics(pgp->pwnd,gp);
	}
	pseries=gp.pfirst;
	while (pseries) {
		Dataseries& series=*pseries;
		plotdata(pgp->pwnd,gp,series);
		pseries=pseries->pnext;
	}
	
	gp.m_update=false;
}


void CGraph::save(FILE*& out,char* name) {

	// Writes data for specified graph (Gp object) to specified i/o stream
	// (normally an ascii text file)

	const int MAXSERIES=100;
	Dataseries* pseries[MAXSERIES];
	Data* pdata[MAXSERIES];
	Gp& gp=getgp(name);
	int s,nseries=0;
	bool matchedx=false;
	Dataseries* p;
	pseries[nseries]=gp.pfirst;
	while (pseries[nseries]) {
		pdata[nseries]=pseries[nseries]->pfirst;
		nseries++;
		if (nseries>=MAXSERIES) {
			fprintf(out,"Too many series for save\n");	
			return;
		}
		pseries[nseries]=pseries[nseries-1]->pnext;
	}
	if (!nseries) {
		fprintf(out,"No data\n");	
		return;
	}

	if (nseries>1) {
	
		// Check if X values match in all series

		bool oksofar=true;
		bool finished=false;

		while (oksofar && !finished) {
			if (pdata[0]) {
				for (s=1;s<nseries;s++) {
					if (pdata[s]) {
						if (pdata[s]->x!=pdata[0]->x) oksofar=false;
						pdata[s]=pdata[s]->pnext;
					}
					else oksofar=false;
				}
				pdata[0]=pdata[0]->pnext;
			}
			else {
				for (s=1;s<nseries;s++)
					if (pdata[s]) oksofar=false;
				finished=true;
			}
		}
		if (finished && oksofar) matchedx=true;
	}

	if (matchedx) {
		fprintf(out,"\t%s\nx-val",name);
		for (s=0;s<nseries;s++) {
			fprintf(out,"\t%s",pseries[s]->name);
			pdata[s]=pseries[s]->pfirst;
		}
		while (pdata[0]) {
			fprintf(out,"\n%g",pdata[0]->x);
			for (s=0;s<nseries;s++) {
				fprintf(out,"\t%g",pdata[s]->y);
				pdata[s]=pdata[s]->pnext;
			}
		}
		fprintf(out,"\n");
	}
	else {
		for (s=0;s<nseries;s++) {
			if (s>0) fprintf(out,"\n");
			fprintf(out,"x-val\t%s-%s\n",name,pseries[s]->name);
			pdata[s]=pseries[s]->pfirst;
			while (pdata[s]) {
				fprintf(out,"%g\t%g\n",pdata[s]->x,pdata[s]->y);
				pdata[s]=pdata[s]->pnext;
			}
		}
	}
}

void CGraph::cleanup() {

	// Removes all Gp and series from this CGraph object

	Gp* pgp,*pnextgp;
	Dataseries* pseries,*pnextseries;
	Data* pdata,*pnextdata;

	if (nwindow) {
	
		pgp=firstwindow;
		while (pgp) {
			pnextgp=pgp->pnext;
			pseries=pgp->pfirst;
			while (pseries) {
				pnextseries=pseries->pnext;
				pdata=pseries->pfirst;
				while (pdata) {
					pnextdata=pdata->pnext;
					delete pdata;
					pdata=pnextdata;
				}
				delete pseries;
				pseries=pnextseries;
			}
			delete pgp;
			pgp=pnextgp;
		}
		nwindow=0;
		firstwindow=NULL;
	}
}

void CGraph::cleargp(Gp* pgpx) {

	// Removes specified Gp and series from this CGraph object

	Gp* pgp,*pnextgp,*plastgp;
	Dataseries* pseries,*pnextseries;
	Data* pdata,*pnextdata;

	if (nwindow) {

		pgp=firstwindow;
		plastgp=NULL;
		while (pgp) {

			if (pgp==pgpx) {
				pnextgp=pgpx->pnext;
				pseries=pgpx->pfirst;
				pgpx->ndata=0;
				while (pseries) {
					pnextseries=pseries->pnext;
					pdata=pseries->pfirst;
					pseries->pfirst=NULL;
					while (pdata) {
						pnextdata=pdata->pnext;
						delete pdata;
						pdata=pnextdata;
					}
					delete pseries;
					pseries=pnextseries;
				}
				delete pgpx;
				if (plastgp)
					plastgp->pnext=pnextgp;
				else
					firstwindow=pnextgp;
				nwindow--;
				return;
			}
			else {
				plastgp=pgp;
				pgp=pgp->pnext;
			}
		}
	}
}

void CGraph::reset() {

	// Removes data from all series and resets Gp objects, but does not delete
	// Graphs. Series are deleted.
	// Applies only to non-tagged windows

	Gp* pgp;
	Dataseries* pseries,*pnextseries;
	Data* pdata,*pnextdata;

	if (nwindow) {
		pgp=firstwindow;
		pgp->ndata=0;
		while (pgp) {
			if (!(pgp->tag)) {
				pseries=pgp->pfirst;
				pgp->pfirst=NULL;
				pgp->ndata=0;
				pgp->nseries=0;
				while (pseries) {
					pnextseries=pseries->pnext;
					pdata=pseries->pfirst;
					pseries->pfirst=NULL;
					while (pdata) {
						pnextdata=pdata->pnext;
						delete pdata;
						pdata=pnextdata;
					}
					delete pseries;
					pseries=pnextseries;
				}
			}

			pgp->m_update=true;
			pgp->pwnd->update_window();

			pgp=pgp->pnext;
		}
	}
}

void CGraph::resetgp(Gp* pgp) {

	// Removes data from a particular Gp object (corresponding to a particular graph window)

	Dataseries* pseries;
	Data* pdata,*pnextdata;

	if (pgp) {
		pseries=pgp->pfirst;
		pgp->ndata=0;
		while (pseries) {
			pdata=pseries->pfirst;
			pseries->pfirst=NULL;
			while (pdata) {
				pnextdata=pdata->pnext;
				delete pdata;
				pdata=pnextdata;
			}
			pseries=pseries->pnext;
		}

		pgp->m_update=true;
		pgp->pwnd->update_window();
	}
}

void CGraph::resetwindow(char* wname) {

	// Removes data from the specified graph window

	if (nwindow) {
		xtring n1,n2;
		n1=wname;
		n1=n1.upper();
		Gp* pgp=firstwindow;
		n2=pgp->name;
		n2=n2.upper();
		if (n1==n2) {
			resetgp(pgp);	
			return;
		}
		while (pgp->pnext) {
			pgp=pgp->pnext;
			n1=wname;
			n1=n1.upper();
			n2=pgp->name;
			n2=n2.upper();
			if (n1==n2) {
				resetgp(pgp);
				return;
			}
		}
	}
	return; // not found
}


///////////////////////////////////////////////////////////////////////////////////////
// FILE SCOPE GLOBAL VARIABLES

FILE* logfile;
pthread_mutex_t mutex;
	// mutex for coordination between model thread and windows manager thread
XClient xclient;
	// X Windows client object
CGraph graph;
	// Graphics object


///////////////////////////////////////////////////////////////////////////////////////
// X WINDOWS MANAGER THREAD

void* listen(void*) {

	// Main thread function for window manager
	// Redraws windows on event messages from X server

	XEvent event;
	bool connected;
	bool redraw;
	Window eventwindow;
	int command;

	pthread_mutex_lock(&mutex);
	connected=xclient.m_connected;
	pthread_mutex_unlock(&mutex);

	while (connected) {

		bool found;

		XNextEvent(xclient.m_display,&event);
		if (event.type==Expose || event.type==ConfigureNotify || event.type==ButtonPress) {
			// Update display of window

			XWindow* pthis;

			// Lock mutex for seek and display
			pthread_mutex_lock(&mutex);

			if (event.type==ConfigureNotify) {
				eventwindow=event.xconfigure.window;
			}
			else if (event.type==Expose) {
				eventwindow=event.xexpose.window;
			}
			else {
				eventwindow=event.xbutton.window;
			}

			pthis=xclient.pfirst;
			found=false;
			while (pthis && !found) {
				if (pthis->m_window==eventwindow) {
					if (!pthis->m_shown) {

						pthis->m_shown=true;
						pthis->set_window_name();

						pthis->m_menu.addstring("Freeze");

						pthis->m_menu.addstring("Clear");
						xtring x1="Save data to ";
						x1+=pthis->m_filename+".txt";

						pthis->m_menu.addstring((char*)x1);

					}
					if (pthis->m_frozen) pthis->m_frozen=false;
					if (event.type==Expose) {
						if (event.xexpose.send_event) {
							redraw=false;
								// Expose events forced by CGraph don't redraw all
						}
						else {
							redraw=true;
								// (system expose events redraw)
						}
					}
					else if (event.type==ConfigureNotify) {
						pthis->resize_window(event.xconfigure.width,event.xconfigure.height);
						pthis->m_pgp->configure(event.xconfigure.width,
							event.xconfigure.height,0);
						pthis->m_menu.showmenu(false,0,0);
						xclient.save_window_specs(WSFILE);
						redraw=true;
					}
					else if (event.type==ButtonPress) {
						if (pthis->m_menu.m_visible) {
							if (event.xbutton.button==Button1) { // chose from menu?
								command=pthis->m_menu.getmenuitem(event.xbutton.x,
									event.xbutton.y);
								pthis->m_menu.showmenu(false,0,0);
								if (command==0) { // Freeze
									pthis->m_frozen=true;
								}
								else if (command==1) {
									graph.resetgp(pthis->m_pgp);	
								}
								else if (command==2) {
									xtring filename=pthis->m_filename+".txt";
									FILE* out=fopen(filename,"wt");
									if (out) {
										graph.save(out,pthis->m_name);
										fclose(out);
										dprintf("Saved %s data to %s\n",
											(char*)pthis->m_name,(char*)filename);
									}
									else {
										dprintf("Could not open %s for output\n",
											(char*)filename);
									}
								}
							}
							else { // move menu
								pthis->m_menu.showmenu(false,0,0);
								pthis->m_menu.showmenu(true,event.xbutton.x,event.xbutton.y);
							}
						}
						else if (event.xbutton.button!=Button1) { // show menu
							pthis->m_menu.showmenu(true,event.xbutton.x,event.xbutton.y);
						}
						redraw=true;
					}

					found=true;

					pthis->draw_window(&graph,redraw);

				}
				else pthis=pthis->pnext;
			}

			connected=xclient.m_connected;
			
			// Unlock mutex
			pthread_mutex_unlock(&mutex);

		}
	}

	return NULL;
};

/*void formatf(xtring& output,char* format,va_list& v) {


	const int MINBUF=100;
	char* pchar;
	xtring sect;
	xtring fspec;
	xtring wspec;
	xtring buffer(MINBUF);
	char* pfmt=(char*)format;
	int i;
	char waitspec=0,readspec=0,readwidth,havewidth;

	output="";
	i=0;
	while (*pfmt) {
		if (waitspec) {
			if (*pfmt>='0' && *pfmt<='9' || *pfmt=='-' || *pfmt=='+' || *pfmt==' ' ||
				*pfmt=='#' || *pfmt=='h' || *pfmt=='l' || *pfmt=='I' || *pfmt=='L' ||
				*pfmt=='c' || *pfmt=='C' || *pfmt=='d' || *pfmt=='i' || *pfmt=='o' ||
				*pfmt=='u' || *pfmt=='x' || *pfmt=='X' || *pfmt=='e' || *pfmt=='E' ||
				*pfmt=='f' || *pfmt=='g' || *pfmt=='G' || *pfmt=='n' || *pfmt=='p' ||
				*pfmt=='s' || *pfmt=='S' || *pfmt=='.') {
				readspec=1;
				havewidth=0;
				readwidth=0;
				wspec="0";
				fspec="%";
				pfmt--;
			}
			else output+=*pfmt;
			waitspec=0;
		}
		else if (readspec) {
			if (*pfmt!='*') fspec+=*pfmt;
			if (readwidth) {
				if (*pfmt>='0' && *pfmt<='9') wspec+=*pfmt;
				else {
					readwidth=0;
					havewidth=1;
				}
			}
			else if (*pfmt>='0' && *pfmt<='9' && !havewidth) {
				wspec=*pfmt;
				readwidth=1;
			}
			else if (*pfmt=='*' && !havewidth) {
				wspec.reserve(MINBUF);
				sprintf(wspec,"%d",va_arg(v,int));
				fspec+=wspec;
				havewidth=1;
			}
			else if (*pfmt!='-' && *pfmt!='+' && *pfmt!=' ' && *pfmt!='#') {
				havewidth=1;
			}
			if (*pfmt=='c' || *pfmt=='C' || *pfmt=='d' || *pfmt=='i' || *pfmt=='o' ||
				*pfmt=='u' || *pfmt=='x' || *pfmt=='X') {
				buffer.reserve(MINBUF+wspec.num());
				sprintf(buffer,fspec,va_arg(v,int));
				output+=buffer;
				readspec=0;
			}
			else if (*pfmt=='e' || *pfmt=='E' || *pfmt=='f' || *pfmt=='g' ||
				*pfmt=='G') {
				buffer.reserve(MINBUF+wspec.num());
				sprintf(buffer,fspec,va_arg(v,double));
				output+=buffer;
				readspec=0;
			}
			else if (*pfmt=='p') {
				buffer.reserve(MINBUF+wspec.num());
				sprintf(buffer,fspec,va_arg(v,void*));
				output+=buffer;
				readspec=0;
			}
			else if (*pfmt=='s' || *pfmt=='S') {
				pchar=va_arg(v,char*);
				buffer.reserve(strlen(pchar)+wspec.num());
				sprintf(buffer,fspec,pchar);
				output+=buffer;
				readspec=0;
			}
		}
		else if (*pfmt=='%') waitspec=1;
		else {
			sect[i++]=*pfmt;
			if (*(pfmt+1)=='%' || !*(pfmt+1)) {
				sect[i]='\0';
				output+=sect;
				i=0;
			}
		}
		pfmt++;
	}
}*/

void fail(xtring format,...) {

        // printf-style function accessible throughout the model code.
        // Sends text to stdio (screen) and log file, then terminates program

        va_list v;
        va_start(v,format);

        xtring output;
        formatf(output,format,v);

        // Produce output
        // (comment out one or both of these if output to the screen and/or log file
        // is not required)

        printf("Fail called\n");
        fprintf(stdout,"%s\n",(char*)output);
        fprintf(logfile,"%s\n",(char*)output);

        exit(99);
}


void dprintf(xtring format,...) {

        // printf-style function accessible throughout the model code.
        // Sends text to stdio (screen) and log file.

        va_list v;
        va_start(v,format);

        xtring output;
        formatf(output,format,v);

        // Produce output
        // (comment out one or both of these if output to the screen and/or log file
        // is not required)

        fprintf(stdout,"%s",(char*)output);
        fprintf(logfile,"%s",(char*)output);
        fflush(logfile);
}


void plot(xtring window_name,xtring series_name,double x,double y) {

	// Adds data point (x,y) to series 'series_name' of line graph 'window_name'. If
	// the series and/or line graph do not yet exist, they are created

	if (ifgra) {
		pthread_mutex_lock(&mutex);
		if (xclient.m_connected) {
			graph.createwindow(&xclient,(char*)window_name,width,height,use_saved_window_dims);
			graph.newseries((char*)window_name,(char*)series_name);
			graph.data((char*)window_name,(char*)series_name,x,y);
		}
		pthread_mutex_unlock(&mutex);
	}
}


void resetwindow(xtring window_name) {

	// 'Forgets' series and data for line graph 'window_name' created using function
	// plot (above)

	if (ifgra) {
		pthread_mutex_lock(&mutex);
		if (xclient.m_connected) graph.resetwindow((char*)window_name);
		pthread_mutex_unlock(&mutex);
	}
}


void clear_all_graphs() {

	// 'Forgets' series and data for all currently-defined line graphs created using
	// function plot (above)

	if (ifgra) {
		pthread_mutex_lock(&mutex);
		if (xclient.m_connected) graph.reset();
		pthread_mutex_unlock(&mutex);
	}
}


bool abort_request_received() {

	// Can't do anything here
	
	return false;
}


///////////////////////////////////////////////////////////////////////////////////////
// HANDLING OF COMMAND-LINE ARGUMENTS

void parseargs(int argc,char* argv[],char* margv[2]) {

	int a,ival;
	xtring word,param;

	int insarg=-1;

	bool failure=false;
	for (a=1;a<argc && !failure;a++) {
		word=argv[a];
		if (word[0]=='-') {
			if (word=="-nogra") {
				ifgra=false;
			}
			else if (word=="-wait") {
				ifwait=true;
			}
			else if (word=="-width") {
				if (a==argc-1) {
					fprintf(stdout,"main: no value specified for -width\n");
					failure=true;
				}
				else {
					param=argv[a+1];
					if (param.isnum()) {
						ival=param.num();
						if (ival<1 || ival>1600) {
							fprintf(stdout,"main: out-of-range value specified for -width\n");
							failure=true;
						}
						else {
							width=ival;
							use_saved_window_dims=false;
							a++;
						}
						
					}
					else {
						fprintf(stdout,"main: qualifier to -width must be a number\n");
						failure=true;
					}
				}
			}
			else if (word=="-height") {
				if (a==argc-1) {
					fprintf(stdout,"main: no value specified for -height\n");
					failure=true;
				}
				else {
					param=argv[a+1];
					if (param.isnum()) {
						ival=param.num();
						if (ival<1 || ival>1200) {
							fprintf(stdout,"main: out-of-range value specified for -height\n");
							failure=true;
						}
						else {
							height=ival;
							use_saved_window_dims=false;
							a++;
						}
						
					}
					else {
						fprintf(stdout,"main: qualifier to -height must be a number\n");
						failure=true;
					}
				}
			}
			else if (word=="-display") {
				if (a==argc-1) {
					fprintf(stdout,"main: nothing specified for -display\n");
					failure=true;
				}
				else {
					display=argv[a+1];
					a++;
				}
			}
			else if (word=="-mono") {
				ifmono=true;
			}
			else if (word=="-log") {
				if (a==argc-1) {
					fprintf(stdout,"main: no filename specified for -log\n");
					failure=true;
				}
				else {
					file_log=argv[a+1];
					a++;
				}
			}
			else if (word=="-nolog") {
				iflog=false;
			}
			else if (word=="-help") {
				insarg=a;
			}
			else {
				fprintf(stdout,"main: invalid option\n");
				failure=true;
			}
		}
		else {
			if (insarg!=-1) {
				fprintf(stdout,"main: more than one instruction file specified\n");
				failure=true;
			}
			else
				insarg=a;
		}
	}

	// Don't wait if there are no graphics anyway
	if (!ifgra) ifwait=false;

	if (insarg==-1) {
		fprintf(stdout,"main: instruction file not specified\n");
		failure=true;
	}

	if (failure) {
		fprintf(stdout,"Usage:   %s [options] <instruction-file>\n",argv[0]);
		fprintf(stdout,"Options: -log <file>       specifies log file name\n");
		fprintf(stdout,"         -nolog            no log file\n");
		fprintf(stdout,"         -display <dest>   address or ip-number to X windows server\n");
		fprintf(stdout,"         -mono             monochrome graphical output\n");
		fprintf(stdout,"         -width <pixels>   graph window width\n");
		fprintf(stdout,"         -height <pixels>  graph window height\n");
		fprintf(stdout,"         -nogra            suppress graphical output\n");
		fprintf(stdout,"         -wait             wait for keyboard at end of run (to hold graphics)\n");
		fprintf(stdout,"         -help             list ins file parameters\n");
		exit(99);
	}

	margv[0]=new char[strlen(argv[0])+1];
	if (margv[0]) strcpy(margv[0],argv[0]);
	else {
		fprintf(stdout,"main: out of memory\n");
		exit(99);
	}
	margv[1]=new char[strlen(argv[insarg])+1];
	if (margv[1]) strcpy(margv[1],argv[insarg]);
	else {
		fprintf(stdout,"main: out of memory\n");
		exit(99);
	}
}


///////////////////////////////////////////////////////////////////////////////////////
// MAIN
// This is the function called when the executable is run

int main(int argc,char* argv[]) {

	char* margv[2];

	// Parse and implement command line arguments
	parseargs(argc,argv,margv);

	// Open log file if possible
	// (comment out if log file output not desired - you will have to comment
	// out the corresponding fprintf's in functions dprintf and fail above also)

	if (iflog) {
		logfile=fopen(file_log,"wt");
		if (!logfile) {
			fprintf(stdout,"main: could not open log file %s for output",
				(char*)file_log);
			exit(99);
		}
	}

	if (ifgra) {

		// Open display

		xclient.initialise((char*)display);
		if (!xclient.m_connected) {
			dprintf("main: could not establish X connection");
			if (display!="")
				dprintf(" to %s",(char*)display);
			dprintf("\n");
		}

		// Initialise mutex (lock variable for coordination between model thread
		// and windows manager thread)

		pthread_t u;
		pthread_attr_t attr;
		pthread_attr_init(&attr);
		pthread_mutexattr_t mattr;
		pthread_mutexattr_init(&mattr);
		pthread_mutex_init(&mutex,&mattr);
		pthread_mutex_unlock(&mutex);
		
		// Start the X windows manager
		pthread_create(&u,&attr,listen,NULL);
	}

	// Call the framework
	framework(2,margv);

	// Cleanup

	delete margv[0];
	delete margv[1];
	graph.cleanup();
	
	// Say goodbye

        if (ifwait) {
		printf("Press [Enter] ");
        	readfor(stdin,"");
	}

	dprintf("\nFinished\n");
	if (logfile && iflog) fclose(logfile);

	return 0;
}
