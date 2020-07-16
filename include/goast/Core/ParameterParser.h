// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef __PARAMSPARSER_H
#define __PARAMSPARSER_H

#include "Auxiliary.h"

//
class Variable {
    
protected:    
  std::string varStr;  
    
public:
  enum { VAR_DOUBLE = 1,
         VAR_INT    = 2,
         VAR_STRING = 3 };

  Variable ( const char *VarStr ) {
    varStr = VarStr;
  }

  Variable() {
    varStr = "";
  }

  int type() const {
    return type ( getVarStr() );
  }

  bool operator== ( const Variable &Var ) const {
    return !strcmp ( Var.getVarStr(), getVarStr() );
  }

  Variable &operator= ( const Variable &Var ) {
    varStr = Var.varStr;
    return *this;
  }

  void dump ( std::ostream &Off = std::cerr ) const {
    Off << varStr;
  }

  double getDouble() const {
    return atof ( getVarStr() );
  }

  int getInt() const {
    return std::atoi ( getVarStr() );
  }

  const char *getVarStr() const {
    return varStr.c_str();
  }

  void setVarStr ( const char *VarStr ) {
    varStr = VarStr;
  }  
  
protected:
  int type ( const char *VarStr ) const {
    char *err;
    std::strtol ( VarStr, &err, 0 );
    if ( *err == '\0' ) return VAR_INT;
    std::strtod ( VarStr, &err );
    if ( *err == '\0' ) return VAR_DOUBLE;
    return VAR_STRING;
  }

};

//
class VariableField {

protected:
  int _numDim;
  int *_dimSizes;
  std::string _nameOfVarField;
  std::vector<Variable*> _field;    
    
public:
  VariableField ( const char *Name ) : _numDim ( 0 ), _dimSizes ( NULL ) {
    _nameOfVarField =  Name;
  }

  VariableField ( int Num ) : _dimSizes ( NULL ) {
    setNumDim ( Num );
  }

  ~VariableField() {
    if ( _dimSizes )
      delete[] _dimSizes;
    std::vector<Variable*>::const_iterator it;
    for ( it = _field.begin(); it != _field.end(); ++it ) {
      delete *it;
    }
  }

  void setNumDim ( int Num ) {
    _numDim = Num;
    if ( _dimSizes )
      delete[] _dimSizes;
    _dimSizes = new int[ Num ];
  }

  int getNumDim() const {
    return _numDim;
  }

  void setDimSize ( int Num, int Size ) {
    if ( Num < _numDim && Num >= 0 ) {
      _dimSizes[ Num ] = Size;
    } else {
      throw BasicException ( "ERROR in VariableField::setDimSize: Num out of range." );
    }
  }

  int getDimSize ( int Num ) const {
    if ( Num < _numDim && Num >= 0 ) {
      return _dimSizes[ Num ];
    } else {
      throw BasicException ( "ERROR in VariableField::getDimSize: Num out of range." );
      return -1;
    }
  }

  void append ( const char *Var ) {
    _field.push_back ( new Variable ( Var ) );
  }

  void dump ( std::ostream &Off = std::cerr ) const {
    Off << "VarField: " << _nameOfVarField << std::endl;
    Off << "numDim = " << _numDim << std::endl;
    std::vector<Variable*>::const_iterator it;
    for ( it = _field.begin(); it != _field.end(); ++it ) {
      ( *it )->dump ( Off );
      Off << ", ";
    }
    Off << std::endl;
    if ( _numDim > 0 ) {
      for ( int i = 0; i < _numDim; i++ )
        Off << "dimSizes[ " << i << " ]= " << _dimSizes[ i ] << std::endl;
    }
  }

  int isSingleField() const {
    return ( _numDim == 0 );
  }

  void getVariable ( int I, Variable &Var ) const {
    if ( I < static_cast<int> ( _field.size() ) ) {
      Var = *_field[ I ];
    } else {
      std::cerr << "index out of bounds\n";
    }
  }

  void getVariable ( int I1, int I2, Variable &Var ) const {
    if ( I1 < _dimSizes[ 0 ] && I2 < _dimSizes[ 1 ] ) {
      Var = *_field[ I1 * _dimSizes[ 1 ] + I2 ];
    } else {
      std::cerr << "index out of bounds\n";
    }
  }

  void getVariable ( Variable &Var ) const {
    if ( _field.size() == 1 ) {
      Var = *_field[ 0 ];
    } else {
      std::cerr << "ERROR in VariableField::getVariable: size = " << _field.size() << std::endl;
    }
  }

  std::vector<Variable*> &getFieldReference ( ) {
    return _field;
  }

  const char *getName() const {
    return _nameOfVarField.c_str();
  }

};


/**
 * \brief A class for reading parameters from configuration files.
 *
 * A parameter file may contain the following type of lines:
 * - parameter line, i.e. PARAMETERNAME VALUE
 * - comment lines, i.e. lines starting with #
 * - white space lines, i.e. lines containing only white spaces
 * - blank lines, i.e. lines containing only EOL
 */
class ParameterParser {
    
    
  std::vector<VariableField*> varFields;

  std::ifstream _inStream;
  bool _echo;
  std::ostream  &_outStream;
  std::string   _outString;
  std::string   _parameterFileName;
  
public:
  //! Default constructor. Since this doesn't parse any parameter file,
  //! you can't do much with an instance created by this constructor.
  ParameterParser ( ) : _echo ( false ), _outStream ( std::cout ) { }

  /*! Instantiate a parser for a certain file. */
  explicit ParameterParser ( const std::string & ParFilename ) : _echo ( false ), _outStream ( std::cout ) {
    initialize ( ParFilename );
  }

  ~ParameterParser() {
    std::vector<VariableField*>::iterator it;
    for ( it = varFields.begin(); it != varFields.end(); ++it )
      delete *it;
  }

  int getDimSize ( const char *VarName, int I = 0 ) const {
    VariableField* var = findFirstVariableField ( VarName );
    return var->getDimSize ( I );
  }

  int getNumDim ( const char *VarName ) const {
    VariableField* var = findFirstVariableField ( VarName );
    return var->getNumDim ();
  }


  double getDouble ( const char *VarName ) const {
    std::vector<VariableField*>::const_iterator it;
    for ( it = varFields.begin(); it != varFields.end(); ++it ) {
      if ( ( strcmp ( ( *it )->getName(), VarName ) == 0 ) && ( *it )->isSingleField() ) {
        Variable Var;
        ( *it )->getVariable ( Var );
        if ( Var.type() == Variable::VAR_DOUBLE || Var.type() == Variable::VAR_INT ) {
          if ( _echo ) _outStream << VarName << " = " << Var.getDouble() << std::endl;
          return Var.getDouble();
        }
      }
    }

    throw BasicException ( "ParameterParser::getDouble(): No match found for double ", VarName );
    return 0.0;
  }
  
  double getDouble ( const char *VarName, int I ) const {
    std::vector<VariableField*>::const_iterator it;
    for ( it = varFields.begin(); it != varFields.end(); ++it ) {
      if ( ( strcmp ( ( *it )->getName(), VarName ) == 0 ) && ( *it )->getNumDim() == 1 ) {
        Variable Var;
        ( *it )->getVariable ( I, Var );
        if ( Var.type() == Variable::VAR_DOUBLE || Var.type() == Variable::VAR_INT ) {
          if ( _echo ) _outStream << VarName << " = " << Var.getDouble() << std::endl;
          return Var.getDouble();
        }
      }
    }
    throw BasicException ( "ParameterParser::getDouble(): No match found for double ", VarName  );
    return 0.0;
  }
  
  double getDoubleOrDefault ( const char *VarName, double Default ) const {
    return hasVariable ( VarName ) ? getDouble ( VarName ) : Default;
  }

  int getInt ( const char *VarName ) const {
    std::vector<VariableField*>::const_iterator it;
    for ( it = varFields.begin(); it != varFields.end(); ++it ) {
      if ( ( strcmp ( ( *it )->getName(), VarName ) == 0 ) && ( *it )->isSingleField() ) {
        Variable Var;
        ( *it )->getVariable ( Var );
        if ( Var.type() == Variable::VAR_INT ) {
          if ( _echo ) _outStream << VarName << " = " << Var.getInt() << std::endl;
          return Var.getInt();
        }
      }
    }
    throw BasicException ( "ParameterParser::getInt(): No match found for integer ", VarName  );
    return 0;
  }
  
  int getInt ( const char *VarName, int I ) const {
    std::vector<VariableField*>::const_iterator it;
    for ( it = varFields.begin(); it != varFields.end(); ++it ) {
      if ( ( strcmp ( ( *it )->getName(), VarName ) == 0 ) && ( *it )->getNumDim() == 1 ) {
        Variable Var;
        ( *it )->getVariable ( I, Var );
        if ( Var.type() == Variable::VAR_INT ) {
          if ( _echo ) _outStream << VarName << " = " << Var.getInt() << std::endl;
          return Var.getInt();
        }
      }
    }
    throw BasicException ( "ParameterParser::getInt(): No match found for integer ",  VarName );
    return 0;
  }
  
  int getIntOrDefault ( const char *VarName, int Default ) const {
    return hasVariable ( VarName ) ? getInt ( VarName ) : Default;
  }  
  
  void getString ( const char *VarName, char *DestStr ) const {
    std::vector<VariableField*>::const_iterator it;
    for ( it = varFields.begin(); it != varFields.end(); ++it ) {
      if ( ( strcmp ( ( *it )->getName(), VarName ) == 0 ) && ( *it )->isSingleField() ) {
        Variable Var;
        ( *it )->getVariable ( Var );
        strcpy ( DestStr, Var.getVarStr() );
        if ( _echo ) _outStream << VarName << " = " << DestStr << std::endl;
        return;
      }
    }
    throw BasicException ( "ParameterParser::getString(): No match found for ", VarName );
  }

  void getString ( const char *VarName, char *DestStr, int I ) const {
    std::vector<VariableField*>::const_iterator it;
    for ( it = varFields.begin(); it != varFields.end(); ++it ) {
      if ( ( strcmp ( ( *it )->getName(), VarName ) == 0 ) && ( *it )->getNumDim() == 1 ) {
        Variable Var;
        ( *it )->getVariable ( I, Var );
        strcpy ( DestStr, Var.getVarStr() );
        if ( _echo ) _outStream << VarName << " = " << DestStr << std::endl;
        return;
      }
    }
    throw BasicException ( "ParameterParser::getString(): No match found for ", VarName );
  }

  std::string getString ( const char *VarName ) const {
    char tempCString[4096];
    getString ( VarName, tempCString );
    return tempCString;
  }

  std::string getString ( const char *VarName, const int I ) const {
    char tempCString[4096];
    getString ( VarName, tempCString, I );
    return tempCString;
  }

  std::string getStringOrDefault ( const char *VarName, std::string Default ) const {
    return hasVariable ( VarName ) ? getString ( VarName ) : Default;
  }

  //! check whether parameter file has a parameter
  bool   hasVariable ( const char* VarName ) const {
    std::vector<VariableField*>::const_iterator it;
    for ( it = varFields.begin(); it != varFields.end(); ++it )
      if ( strcmp ( ( *it )->getName(), VarName ) == 0 )
        return true;

    // variable name not found:
    return false;
  }


  //! check whether parameter file has a variable, print
  //! error message to cerr if not
  bool   checkVariable  ( const char * VarName ) const {
    if ( hasVariable ( VarName ) )
      return true;

    std::cerr << "ParameterParser::checkVariable(): Parameter file \"" << _parameterFileName << "\" is supposed "
          "to contain field \"" << VarName << "\"." << std::endl;
    return false;
  }

  //! Returns true if the parameter file has the variable \arg VarName and its int value is 1.
  //! Otherwise, false is returned.
  bool   checkAndGetBool ( const char *VarName ) const {
    return ( hasVariable ( VarName ) && ( getInt ( VarName ) == 1 ) );
  }

  //! Returns true/false if the parameter file has the variable \arg VarName and its int value is 1/0.
  //! Throws an exception if either the variable doesn't exist or has a value different from 0 and 1.
  bool   getBool ( const char *VarName ) const {
    const int value = getInt ( VarName );
    if ( value == 1 )
      return true;
    else if ( value == 0 )
      return false;

    throw BasicException ( "ParameterParser::getBool(): invalid value for bool. Should be either 0 or 1." );
    return false;
  }
  
  //! Returns true/false if the parameter file has the variable \arg VarName and its int value is 1/0
  //! Otherwise, returns \arg Default
  bool getBoolOrDefault ( const char *VarName, bool Default ) const {
    return hasVariable ( VarName ) ? getBool ( VarName ) : Default;
  }
  
  void   getIntVec ( const char *VarName, std::vector<int> &vec ) const {
    const int size = this->getDimSize ( VarName );
    vec.resize ( size );
    for ( int i = 0; i < size; ++i )
      vec[i] = this->getInt ( VarName, i );
  }
  
  template<typename RealType>
  void getRealVec ( const char *VarName, std::vector<RealType> &vec ) const {
    const int size = this->getDimSize ( VarName );
    vec.resize ( size );
    for ( int i = 0; i < size; ++i ) 
      vec[i] = static_cast<RealType>( this->getDouble ( VarName, i ) );
  }

  const std::vector<VariableField*>& getVariableFields() const {
    return varFields;
  }

  // initialize parser
  void initialize ( const std::string & ParFilename ) {
    _inStream.open ( ParFilename.c_str() );
    if ( _inStream.fail() )
      throw BasicException ( "ParameterParser: Can't open file ", ParFilename  );
    parse();
    _inStream.close();
    _parameterFileName = ParFilename;
} 
  
protected:
VariableField* findFirstVariableField ( const char *VarName ) const {
  std::vector<VariableField*>::const_iterator it;
  for ( it = varFields.begin(); it != varFields.end(); ++it ) {
    if ( ( strcmp ( ( *it )->getName(), VarName ) == 0 ) ) {
      return (*it);
    }
  }
  throw BasicException ( "ParameterParser::findFirstVariableField(): No match found for ", VarName );
  return NULL;
}
  
  
  void parse(){
  char name[ 80 ];
  char Var[ 8*4096 ];

  while ( !_inStream.eof() ) {

    ignoreWS();

    while ( _inStream.peek() == '#' || ( !_inStream.eof() && ws ( _inStream.peek() ) ) || _inStream.peek() == '\n' ) {
      if ( ws ( _inStream.peek() ) )
        ignoreWS();
      else
        ignoreLine();
    }

    if ( !_inStream.eof() ) {
      _inStream >> name;

      // Don't allow any variable to be defined more than once.
      if ( hasVariable ( name ) ){
        std::ostringstream errorMessage;
        errorMessage << "ParameterParser::parse(): Variable \"" << name << "\" already defined";
        throw BasicException ( errorMessage.str() );
      }

      VariableField *varField = new VariableField ( name );

      ignoreWS();

      if ( _inStream.eof() || ( _inStream.peek() == '\n' ) ){
        throw BasicException ( "ParameterParser::parse(): Unexpected end of line while parsing variable ", name );
      }
      else if ( _inStream.peek() == '{' )
        readField ( *varField );
      else {
        varField->setNumDim ( 0 );
        // If enclosed by quotation marks, allow the variable value to contain spaces.
        if ( _inStream.peek() == '"' ) {
          _inStream.ignore();
          std::string fullVar;
          while ( !_inStream.eof() && ( _inStream.peek() != '"' ) && ( _inStream.peek() != '\n' ) ) {
            fullVar += _inStream.get();
          }
          if ( _inStream.peek() == '"' )
            _inStream.ignore();
          else{
            throw BasicException ( "ParameterParser::parse(): Unexpected end of line while parsing quotation mark enclosed variable ", name );
          }
          varField->append ( fullVar.c_str() );
        }
        else {
          _inStream >> Var;
          varField->append ( Var );
        }
      }
      varFields.push_back ( varField );
      //varField->dump( );

      // Ignore any trailing white space and make sure that nothing but white space is left in this line.
      ignoreWS();
      if ( !_inStream.eof() && ( _inStream.peek() != '\n' ) ){
        std::ostringstream errorMessage;
        errorMessage << "ParameterParser::parse():Error while parsing variable  \"" << name << "\": Unexpected data after the variable value found";
        throw BasicException ( errorMessage.str() );
      }
    }
  }
}
  
  
void readField ( VariableField &VarField ){
  std::string varStr;
  int cd = 0;
  int maxd = -1;
  int i, p;

  int dimSizes[ 16 ], prevDimSizes[ 16 ];
  for ( i = 0; i < 16; i++ ) {
    dimSizes[ i ] = 0;
    prevDimSizes[ i ] = -1;
  }

  ignoreWS();

  do {
    p = _inStream.get();
    if ( ws ( p ) || p == '}' ) {
      // We wrote something to varStr, give the string to VarField and increase dimSizes[ cd ].
      if ( varStr.empty() == false ) {
        VarField.append ( varStr.c_str() );
        dimSizes[ cd ]++;
      }
      varStr.clear();
      if ( p == '}' ) {
        if ( maxd < 0 ) {
          maxd = cd;
        }
        // After closing a bracket, check if the bracket has the same size as the previous in this depth.
        if ( dimSizes[ cd ] != 0 ) {
          if ( prevDimSizes[ cd ] == -1 )
            prevDimSizes[ cd ] = dimSizes[ cd ];
          else if ( dimSizes[ cd ] != prevDimSizes[ cd ] ) {
            std::ostringstream errorMessage;
            errorMessage << "ParameterParser::readField(): ERROR sizes in depth " << cd << " not constant! prevDim = " << prevDimSizes[ cd ] << ", dim = " << dimSizes[ cd ];
            throw BasicException ( errorMessage.str() );
          }
        }
        cd--;
        dimSizes[ cd ]++;
      }
      ignoreWS();
    } else if ( p == '{' ) {
      if ( maxd > 0 && cd == maxd ) {
        throw BasicException ( "ParameterParser::readField(): exceeding maximum depth\n" );
      }
      cd++;
      dimSizes[ cd ] = 0;
      ignoreWS();
    } else {
      varStr += static_cast<char> ( p );
    }
  } while ( cd > 0 && !_inStream.eof() );

  if ( cd != 0 ) {
    std::cerr << "ParameterParser::readField(): missing closing \'}\'\n";
  }

  VarField.setNumDim ( maxd );
  for ( i = 0; i < maxd; i++ ) {
    VarField.setDimSize ( i, dimSizes[ i + 1 ] );
  }
}

  int ws ( int c ) const {
    return ( ( c == ' ' || c <= 0x17 ) && ( c != '\n' ) );
  }

  void ignoreWS() {
    while ( !_inStream.eof() && ws ( _inStream.peek() ) ) {
      _inStream.ignore();
    }
  }

  void ignoreLine() {
    while ( !_inStream.eof() && _inStream.peek() != '\n' ) {
      _inStream.ignore();
    }
    if ( !_inStream.eof() && _inStream.peek() == '\n' )
      _inStream.ignore();
  }

}; 

#endif