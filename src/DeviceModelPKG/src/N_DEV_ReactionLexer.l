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
class N_DEV_ReactionLexer;
// now we can include this without errors from bison 2.4
#include "N_DEV_ReactionParser.hxx"
// BLEAH!   This is here DUPLICATED from N_DEV_ReactionParser.yxx
// because of STUPID choice in Bison 2.3 to put the post-prologue into the
// .cxx file instead of the .hxx file that Bison 2.1 used to put it in.
#undef yyFlexLexer
  /* CAREFUL watch continuations! */
#define YY_DECL \
  int N_DEV_ReactionLexer::getToken(N_DEV::N_DEV_ReactionParser::semantic_type *lvalp, \
                                       N_DEV::location *llocp, \
                                       map<string,int> &theSpeciesIDs)

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
#define yylex N_DEVlex
//#include <N_DEV_ReactionLexer.h>
#include "location.hh"
#include "position.hh"

%}
%option noyywrap
%option yyclass="N_DEV_ReactionLexer"
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
        return N_DEV::N_DEV_ReactionParser::token::RXN_INT;
        }

[\-]?{DIGIT}*"."{DIGIT}* {
        lvalp->dval=strtod(YYText(),0);
        return N_DEV::N_DEV_ReactionParser::token::RXN_FLOAT;
        }

[\-]?{DIGIT}*"."?{DIGIT}*[eE][\-\+]?{DIGIT}{1,3} {
        lvalp->dval=strtod(YYText(), 0);
        return N_DEV::N_DEV_ReactionParser::token::RXN_FLOAT;
       }

<INITIAL>species {
  BEGIN(speciedef);
  return (N_DEV::N_DEV_ReactionParser::token::SPECIES_BEGIN);
  }

<INITIAL>reactions {
  BEGIN(reactdef);
  return (N_DEV::N_DEV_ReactionParser::token::REACTIONS_BEGIN);
  }

<INITIAL>sources {
  BEGIN(sourcedef);
  return (N_DEV::N_DEV_ReactionParser::token::SOURCES_BEGIN);
  }

<INITIAL>initial_conditions {
  BEGIN(ics_def);
  return (N_DEV::N_DEV_ReactionParser::token::ICS_BEGIN);
  }

end {
  BEGIN(INITIAL);
  return (N_DEV::N_DEV_ReactionParser::token::KEY_END);
  }

electron_capture  { return(N_DEV::N_DEV_ReactionParser::token::KEY_ECAPTURE);}
electron_emission { return(N_DEV::N_DEV_ReactionParser::token::KEY_EEMISSION);}
hole_capture  { return(N_DEV::N_DEV_ReactionParser::token::KEY_HCAPTURE);}
hole_emission { return(N_DEV::N_DEV_ReactionParser::token::KEY_HEMISSION);}
simple   { return(N_DEV::N_DEV_ReactionParser::token::KEY_SIMPLE);}
complex  { return(N_DEV::N_DEV_ReactionParser::token::KEY_COMPLEX);}
decomplex  { return(N_DEV::N_DEV_ReactionParser::token::KEY_DECOMPLEX);}
capture  { return(N_DEV::N_DEV_ReactionParser::token::KEY_CAPTURE);}
emission { return(N_DEV::N_DEV_ReactionParser::token::KEY_EMISSION);}
const { return(N_DEV::N_DEV_ReactionParser::token::KEY_CONST);}

<speciedef>{ID}   {
     lvalp->sval=new std::string(YYText());
    return N_DEV::N_DEV_ReactionParser::token::SPECIE;
     }

<reactdef,sourcedef,ics_def>{ID}   {
  lvalp->sval=new std::string(YYText());
     if (theSpeciesIDs.find(YYText()) != theSpeciesIDs.end())
     { 
       return N_DEV::N_DEV_ReactionParser::token::SPECIE;
     }
     else
     {
       return N_DEV::N_DEV_ReactionParser::token::IDENTIFIER;
     }
}

<sourcedef>"{"[^}]*"}" {
  lvalp->sval=new std::string(YYText()); 
  (lvalp->sval)->erase(0,1);// lop off open brace
  (lvalp->sval)->erase((lvalp->sval)->length()-1); // lop of close brase
  return N_DEV::N_DEV_ReactionParser::token::XYCE_EXPRESSION;
}

[+;,()]   {
     return((YYText())[0]);
    }


"->" {
     return(N_DEV::N_DEV_ReactionParser::token::ARROW);
    }
 

"#"[^\n\r]*     /* eat up one-line comments */

[ \t]+        llocp->step();  /* eat up whitespace */
(\n+|(\r\n)+|\r+)         llocp->lines(YYLeng());

.           std::cout << "Unrecognized character: " <<  YYText() << std::endl;

%%
