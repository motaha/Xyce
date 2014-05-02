/*
 * OFFICIAL USE ONLY
 * May be exempt from public release under the Freedom of Information Act 
 * (5 U.S.C. 552), exemption number and category: 3 Statutory Exemption
 * Department of Energy review required before public release.
 * Name/Org: Scott Hutchinson, 1437    Date: 5/9/07

 * EXPORT CONTROLLED INFORMATION: 
 * Treat this material per Department of State (DOS) International Traffic in 
 * Arms Regulations, 22CFR120-130.  Information contained in this document is 
 * also subject to controls defined by the Department of Defense Directive 
 * 5230.25. 
//
*/
%{
#include <Xyce_config.h>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <map>
// these here because Bison 2.4 no longer puts the pre-prologue into the
// header file... grrrr.
#include <vector>
#include <N_DEV_Specie.h>
#include <N_DEV_ReactionNetwork.h>
namespace Xyce {
namespace Device {
class ReactionLexer;
}}
// now we can include this without errors from bison 2.4
#include "N_DEV_ReactionParser.hxx"
// BLEAH!   This is here DUPLICATED from ReactionParser.yxx
// because of STUPID choice in Bison 2.3 to put the post-prologue into the
// .cxx file instead of the .hxx file that Bison 2.1 used to put it in.
#undef yyFlexLexer
/* CAREFUL watch continuations! */
#undef YY_DECL
#define YY_DECL \
  int Xyce::Device::ReactionLexer::getToken(XyceDevice::ReactionParser::semantic_type *lvalp, \
                                            XyceDevice::location *llocp, \
                                            std::map<std::string, int> &theSpeciesIDs)

  // YECH!  Work around very stupid way that multiple parsers/lexers are 
  // handled.
  // Bison's "%name-prefix" is implemented as a #define yylex "prefix"lex
  // which BREAKS flex's C++ lexer: it contains a method "yylex" in the
  // yyFlexLexer class.  Unless we do this kludge, that method gets renamed
  // with the define as well, and the result is a broken set of classes
#undef yylex
#include <FlexLexer.h>
#include <N_DEV_ReactionLexer.h>
  // undo that kludge.  Note that because of this stupidity, if the 
  // "%name-prefix" is changed, this line needs to be changed, too.
#define yylex XyceDevicelex
//#include <ReactionLexer.h>
#include "location.hh"
#include "position.hh"

%}
%option noyywrap
%option yyclass="Xyce::Device::ReactionLexer"
%s speciedef
%s reactdef
%s sourcedef
%s ics_def

DIGIT [0-9]
ID [a-zA-Z][a-zA-Z0-9_]*

%{
#define YY_USER_ACTION llocp->columns(YYLeng());
%}

%%
%{
  llocp->step();
%}

[\-]?{DIGIT}+ {
        lvalp->ival=strtol(YYText(),0,10);
        return XyceDevice::ReactionParser::token::RXN_INT;
        }

[\-]?{DIGIT}*"."{DIGIT}* {
        lvalp->dval=strtod(YYText(),0);
        return XyceDevice::ReactionParser::token::RXN_FLOAT;
        }

[\-]?{DIGIT}*"."?{DIGIT}*[eE][\-\+]?{DIGIT}{1,3} {
        lvalp->dval=strtod(YYText(), 0);
        return XyceDevice::ReactionParser::token::RXN_FLOAT;
       }

<INITIAL>species {
  BEGIN(speciedef);
  return (XyceDevice::ReactionParser::token::SPECIES_BEGIN);
  }

<INITIAL>reactions {
  BEGIN(reactdef);
  return (XyceDevice::ReactionParser::token::REACTIONS_BEGIN);
  }

<INITIAL>sources {
  BEGIN(sourcedef);
  return (XyceDevice::ReactionParser::token::SOURCES_BEGIN);
  }

<INITIAL>initial_conditions {
  BEGIN(ics_def);
  return (XyceDevice::ReactionParser::token::ICS_BEGIN);
  }

end {
  BEGIN(INITIAL);
  return (XyceDevice::ReactionParser::token::KEY_END);
  }

electron_capture  { return(XyceDevice::ReactionParser::token::KEY_ECAPTURE);}
electron_emission { return(XyceDevice::ReactionParser::token::KEY_EEMISSION);}
hole_capture  { return(XyceDevice::ReactionParser::token::KEY_HCAPTURE);}
hole_emission { return(XyceDevice::ReactionParser::token::KEY_HEMISSION);}
simple   { return(XyceDevice::ReactionParser::token::KEY_SIMPLE);}
complex  { return(XyceDevice::ReactionParser::token::KEY_COMPLEX);}
decomplex  { return(XyceDevice::ReactionParser::token::KEY_DECOMPLEX);}
capture  { return(XyceDevice::ReactionParser::token::KEY_CAPTURE);}
emission { return(XyceDevice::ReactionParser::token::KEY_EMISSION);}
const { return(XyceDevice::ReactionParser::token::KEY_CONST);}

<speciedef>{ID}   {
     lvalp->sval=new std::string(YYText());
    return XyceDevice::ReactionParser::token::SPECIE;
     }

<reactdef,sourcedef,ics_def>{ID}   {
  lvalp->sval=new std::string(YYText());
     if (theSpeciesIDs.find(YYText()) != theSpeciesIDs.end())
     { 
       return XyceDevice::ReactionParser::token::SPECIE;
     }
     else
     {
       return XyceDevice::ReactionParser::token::IDENTIFIER;
     }
}

<sourcedef>"{"[^}]*"}" {
  lvalp->sval=new std::string(YYText()); 
  (lvalp->sval)->erase(0,1);// lop off open brace
  (lvalp->sval)->erase((lvalp->sval)->length()-1); // lop of close brase
  return XyceDevice::ReactionParser::token::XYCE_EXPRESSION;
}

[+;,()]   {
     return((YYText())[0]);
    }


"->" {
     return(XyceDevice::ReactionParser::token::ARROW);
    }
 

"#"[^\n\r]*     /* eat up one-line comments */

[ \t]+        llocp->step();  /* eat up whitespace */
(\n+|(\r\n)+|\r+)         llocp->lines(YYLeng());

.           std::cout << "Unrecognized character: " <<  YYText() << std::endl;

%%
