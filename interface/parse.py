from enum import Enum
import os
import re
import sys
from typing import Dict, List


################################################################################
#                                                                              #
# Resources                                                                    #
#                                                                              #
################################################################################


FCharLenTDefinition = '''#if defined(TINKER_GFORTRAN) && (__GNUC__ <= 7)
// https://gcc.gnu.org/onlinedocs/gfortran/Argument-passing-conventions.html
typedef int tinker_fchar_len_t;
#else
#include <string.h>
typedef size_t tinker_fchar_len_t;
#endif

typedef struct tinker_fchars_st { char* string; tinker_fchar_len_t capacity; } tinker_fchars;
'''

FortranRuntimeHeader = '''void tinkerFortranRuntimeBegin(int, char**);
void tinkerFortranRuntimeEnd();
'''

FortranRuntimeSrc = '''#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef TINKER_GFORTRAN
void _gfortran_set_args(int, char**);
void tinkerFortranRuntimeBegin(int argc, char** argv) { _gfortran_set_args(argc, argv); }
void tinkerFortranRuntimeEnd() {}
#endif

#ifdef TINKER_IFORT
void for_rtl_init_(int*, char**);
int for_rtl_finish_();
extern int for__l_argc;
extern char** for__a_argv;
void tinkerFortranRuntimeBegin(int argc, char** argv) { for__l_argc = argc; for__a_argv = argv; for_rtl_init_(&argc, argv); }
void tinkerFortranRuntimeEnd() { for_rtl_finish_(); }
#endif

#ifdef __cplusplus
}
#endif'''

ShellPart1: str = '''#!/bin/bash
HEADER={ModHeaderName}.h
MODCPP={ModHeaderName}.cpp
DIR=detail

mkdir -p $DIR
rm -f $HEADER
rm -f $DIR/*.h*

##########

cat << ENDOFFILE >> $HEADER
#pragma once

ENDOFFILE
'''

ShellPart2: str = '''
##########

cat << ENDOFFILE >> $HEADER
ENDOFFILE
'''

CMAKE_LIBTINKER = '''cmake_minimum_required (VERSION 3.0)

project (TinkerInterface LANGUAGES NONE)

enable_language (Fortran)
add_library (tinkerObjF OBJECT
{FORT_FILES})

if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/cpp/tinker/routines.h")
    enable_language (CXX)
    add_library (tinkerObjCpp OBJECT cpp/tinker/routines.cpp cpp/tinker/modcpp.cpp)
    if (${CMAKE_Fortran_COMPILER_ID} STREQUAL GNU)
        target_include_directories (tinkerObjCpp PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/tinker/gfortran")
    elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)
        target_include_directories (tinkerObjCpp PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/tinker/ifort")
    else ()
        message (FATAL_ERROR "Must use a GNU or Intel Fortran compiler; Please export FC=valid_fortran_compiler; ${CMAKE_Fortran_COMPILER_ID} is not supported.")
    endif ()
    add_library (tinkerFToCpp STATIC $<TARGET_OBJECTS:tinkerObjF> $<TARGET_OBJECTS:tinkerObjCpp>)
    target_include_directories (tinkerFToCpp PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/cpp")
    target_link_libraries (tinkerFToCpp PUBLIC tinkerObjF tinkerObjCpp)
endif ()
if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/c/tinker/routines.h")
    enable_language (C)
    add_library (tinkerObjC OBJECT c/tinker/routines.c)
    if (${CMAKE_Fortran_COMPILER_ID} STREQUAL GNU)
        target_include_directories (tinkerObjC PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/tinker/gfortran")
    elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)
        target_include_directories (tinkerObjC PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/tinker/ifort")
    else ()
        message (FATAL_ERROR "Must use a GNU or Intel Fortran compiler; Please export FC=valid_fortran_compiler; ${CMAKE_Fortran_COMPILER_ID} is not supported.")
    endif ()
    add_library (tinkerFToC STATIC $<TARGET_OBJECTS:tinkerObjF> $<TARGET_OBJECTS:tinkerObjC>)
    target_include_directories (tinkerFToC PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/c")
    target_link_libraries (tinkerFToC PUBLIC tinkerObjF tinkerObjC)
endif ()'''


################################################################################
#                                                                              #
# Files                                                                        #
#                                                                              #
################################################################################


class FileType(Enum):
    unknown = 0
    module = 1
    subroutine = 22
    program = 3

    t9_extern_ref = 91
    t9_extern_c = 92
    t9_extern_c_with_macro_const = 93
    t9_definition = 94


class RawFile(List[str]):
    @staticmethod
    def __readLines(filename: str) -> List[str]:
        if filename == '':
            return []

        content: List[str] = [line.rstrip().lower() for line in open(filename)]
        # remove ampersand
        totl: int = len(content)
        if totl == 0 or totl == 1:
            return content

        lines: List[str] = []
        for line in content:
            if len(line) > 6 and line[5] in ['&', '$']:
                lines[-1] = lines[-1] + line[6:]
            else:
                lines.append(line)
        return lines

    def read(self, filename: str) -> None:
        self.clear()
        for x in self.__readLines(filename):
            self.append(x)

    def __init__(self, filename: str = '') -> None:
        super().__init__()
        self.read(filename)

    def __str__(self) -> str:
        return '\n'.join(self)


def getDirnameStemExt(filename: str):
    dirname: str = os.path.dirname(filename)
    if dirname == '':
        dirname = '.'
    base: str = os.path.basename(filename)
    stem, ext = os.path.splitext(base)
    return dirname, stem, ext


class FileInfo:
    def __init__(self) -> None:
        self.file_type: FileType = FileType.unknown
        self.filename: str = ''
        self.raw: RawFile = RawFile()


class FileInfoLists:
    def __init__(self, filenames: List[str]) -> None:
        self.unknown_files: List[FileInfo] = []
        self.subroutine_files: List[FileInfo] = []
        self.program_files: List[FileInfo] = []
        self.module_files: List[FileInfo] = []
        for filename in filenames:
            self.__categorize_file(filename)

    def __categorize_file(self, fortran_filename: str) -> None:
        def type_must_be_unknown(t):
            if t != FileType.unknown:
                raise ValueError('Cannot parse file: {}'.format(fortran_filename))

        content: RawFile = RawFile(fortran_filename)
        file_type = FileType.unknown
        for line in content:
            # '      module foo'
            if len(line) > 13:
                if line[0:13] == '      module ':
                    type_must_be_unknown(file_type)
                    file_type = FileType.module
            # '      program foo'
            if len(line) > 14:
                if line[0:14] == '      program ':
                    type_must_be_unknown(file_type)
                    file_type = FileType.program
        if file_type == FileType.unknown:
            file_type = FileType.subroutine

        finfo = FileInfo()
        finfo.filename = fortran_filename
        finfo.raw = content
        finfo.file_type = file_type
        if file_type == FileType.program:
            self.program_files.append(finfo)
        elif file_type == FileType.subroutine:
            self.subroutine_files.append(finfo)
        elif file_type == FileType.module:
            self.module_files.append(finfo)


################################################################################
#                                                                              #
# Parser                                                                       #
#                                                                              #
################################################################################


class Tinker8:
    @staticmethod
    def fortranExternalCallInC(func: str, arg: str) -> str:
        lookupTable = {}

        fname = 'evalue'
        for routine in ['numgrad']:
            lookupTable[fname + '@' + routine] = 'double (*evalue)()'

        fname = 'fvalue'
        for routine in ['simplex']:
            lookupTable[fname + '@' + routine] = 'double (*fvalue)(double*)'

        fname = 'gvalue'
        for routine in ['bsstep', 'diffeq', 'mmid']:
            lookupTable[fname + '@' + routine] = 'void (*gvalue)(double*, double*, double*)'

        fname = 'fgvalue'
        for routine in ['lbfgs', 'ocvm', 'search', 'tncg', 'tnsolve']:
            lookupTable[fname + '@' + routine] = 'double (*fgvalue)(double*, double*)'

        fname = 'optsave'
        for routine in ['lbfgs', 'ocvm', 'tncg']:
            lookupTable[fname + '@' + routine] = 'void (*optsave)(int*, double*, double*)'

        fname = 'rsdvalue'
        for routine in ['square', 'trust']:
            lookupTable[fname + '@' + routine] = 'void (*rsdvalue)(int*, int*, double*, double*)'

        fname = 'lsqwrite'
        for routine in ['square']:
            lookupTable[fname + '@' + routine] = 'void (*lsqwrite)(int*, int*, double*, double*, double*)'

        fname = 'hmatrix'
        for routine in ['tncg']:
            lookupTable[fname + '@' + routine] = 'void (*hmatrix)(char*, double*, double*, int*, int*, int*, double*, tinker_fchar_len_t)'

        k = arg + '@' + func
        if k in lookupTable.keys():
            return lookupTable[k]
        else:
            raise ValueError('Cannot find the type of {}'.format(k))


def static_initialize(**kwargs):
    def decorate(func):
        for k in kwargs:
            setattr(func, k, kwargs[k])
        return func
    return decorate


@static_initialize(counter=0, content={
    'integer': 'int',
    'integer*8': 'unsigned long long',  # used as memory handle
    'logical': 'int',
    'real*8': 'double',
    'real*4': 'float',
    'external': 'fptr',  # function pointer
    'character*(*)': 'char',  # used as the type of function arguments
    'character': 'char',
    'character*1': 'char'})
def KnownTypes() -> Dict[str, str]:
    if KnownTypes.counter == 0:
        KnownTypes.counter = KnownTypes.counter + 1
        for i in range(2, 1025):
            KnownTypes.content['character*{}'.format(i)] = 'char'
    return KnownTypes.content


def CppKeywords() -> List[str]:
    return [
        'alignas', 'alignof', 'and', 'and_eq', 'asm',
        'atomic_cancel', 'atomic_commit', 'atomic_noexcept', 'auto',
        'bitand', 'bitor', 'bool', 'break',
        'case', 'catch', 'char', 'char8_t', 'char16_t', 'char32_t',
        'class', 'compl', 'concept', 'const', 'consteval', 'constexpr',
        'const_cast', 'continue', 'co_await', 'co_return', 'co_yield',
        'decltype', 'default', 'delelte', 'do', 'double', 'double_cast',
        'else', 'enum', 'explicit', 'export', 'extern',
        'false', 'float', 'for', 'friend',
        'goto',
        'if', 'import', 'inline', 'int',
        'long',
        'module', 'mutable',
        'namespace', 'new', 'noexcept', 'not', 'not_eq', 'nullptr',
        'operator', 'or', 'or_eq',
        'private', 'protected', 'public',
        'reflexpr', 'register', 'reinterpret_cast', 'requires', 'return',
        'short', 'signed', 'sizeof', 'static', 'static_assert',
        'static_cast', 'struct', 'switch', 'synchronized',
        'template', 'this', 'thread_local', 'throw',
        'true', 'try', 'typedef', 'typeid', 'typename',
        'union', 'unsigned', 'using',
        'virtual', 'void', 'volatile',
        'wchar_t', 'while',
        'xor', 'xor_eq']


class Variable:
    def __init__(self) -> None:
        self.dimension: List[str] = []
        self.scope: str = ''
        self.symbol: str = ''  # module atoms; subroutine elj; int n;
        self.value: str = ''
        self.type: str = ''  # module, fptr, int, char, etc.
        self.return_type: str = ''
        self.is_const: bool = False
        self.entries: List[Variable] = []  # function args, module variables.
        # module
        self.module_has_symbol: bool = False
        self.module_depends_on: List[str] = []

    def find_entry_by_symbol(self, symbol: str):
        for e in self.entries:
            if e.symbol == symbol:
                return True, e
        return False, Variable()

    def f2c(self, op: FileType):
        MacroTinkerMod: str = 'TINKER_MOD'

        macro_const: bool = False
        is_c_header: bool = False
        if op == FileType.t9_extern_c_with_macro_const:
            macro_const = True
            is_c_header = True
            op = FileType.t9_extern_c

        extern_C_str: str = 'extern "C"'
        if is_c_header:
            extern_C_str = 'extern'

        line: str = '#error ERROR'
        local_symbol = self.symbol
        if local_symbol in CppKeywords():
            local_symbol = local_symbol + '_'
        if self.is_const:
            # sizes.f: integer maxatm
            if op == FileType.t9_extern_ref:
                # const int a = 100;
                line = 'const {} {} = {};'.format(self.type, local_symbol, self.value)
            elif macro_const and op == FileType.t9_extern_c:
                # #define a 100
                line = '#define {}__{} {}'.format(MacroTinkerMod, self.symbol, self.value)
            else:
                line = ''
        elif len(self.dimension):
            if ':' not in self.dimension:
                # angpot.f: character*8 opbtyp
                # extern int (&a)[maxatm]; <-> integer a(maxatm)
                # extern char (&b)[6][5][4][240]; <-> character*240 b(4,5,6)
                if op == FileType.t9_extern_ref:
                    # extern double (&x)
                    line = 'extern {} (&{})'.format(self.type, local_symbol)
                elif op == FileType.t9_extern_c:
                    # extern "C" double TINKER_MOD(atoms, x)
                    line = '{} {} {}({}, {})'.format(extern_C_str, self.type, MacroTinkerMod, self.scope, self.symbol)
                elif op == FileType.t9_definition:
                    # double (&x)
                    line = '{} (&{})'.format(self.type, local_symbol)
                for x in reversed(self.dimension):
                    if macro_const:
                        dimstr: str = ''
                        try:
                            # [42]
                            dimstr = '[{}]'.format(int(x))
                        except ValueError:
                            # [TINKER_MOD__maxatm]
                            dimstr = '[{}__{}]'.format(MacroTinkerMod, x)
                        line = line + dimstr
                    else:
                        line = line + '[{}]'.format(x)
                if op == FileType.t9_definition:
                    line = line + ' = {}({}, {})'.format(MacroTinkerMod, self.scope, self.symbol)
                line = line + ';'
            elif self.dimension[0] == ':':
                # align.f: integer, allocatable :: ifit(:,:)
                # extern float*& c; <-> real*4, allocatable :: c(:,:)
                if op == FileType.t9_extern_ref:
                    line = 'extern {}*& {};'.format(self.type, local_symbol)
                elif op == FileType.t9_extern_c:
                    line = '{} {}* {}({}, {});'.format(extern_C_str, self.type, MacroTinkerMod, self.scope, self.symbol)
                elif op == FileType.t9_definition:
                    line = '{}*& {} = {}({}, {});'.format(self.type, local_symbol, MacroTinkerMod, self.scope, self.symbol)
            else:
                # angpot.f: character*8, allocatable :: angtyp(:)
                # extern char (*&d)[6]; <-> character*6, allocatable :: d(:,:)
                known_dimension: List[str] = [x for x in self.dimension if x != ':']
                if op == FileType.t9_extern_ref:
                    line = 'extern {} (*&{})'.format(self.type, local_symbol)
                elif op == FileType.t9_extern_c:
                    line = '{} {} (*{}({}, {}))'.format(extern_C_str, self.type, MacroTinkerMod, self.scope, self.symbol)
                elif op == FileType.t9_definition:
                    line = '{} (*&{})'.format(self.type, local_symbol)
                for x in reversed(known_dimension):
                    line = line + '[{}]'.format(x)
                if op == FileType.t9_definition:
                    line = line + ' = {}({}, {})'.format(MacroTinkerMod, self.scope, self.symbol)
                line = line + ';'
        else:
            # action.f: integer neb
            if op == FileType.t9_extern_ref:
                # extern int& n;
                line = 'extern {}& {};'.format(self.type, local_symbol)
            elif op == FileType.t9_extern_c:
                # extern "C" int TINKER_MOD(atoms, n);
                line = '{} {} {}({}, {});'.format(extern_C_str, self.type, MacroTinkerMod, self.scope, self.symbol)
            elif op == FileType.t9_definition:
                # int& n = TINKER_MOD(atoms, n);
                line = '{}& {} = {}({}, {});'.format(self.type, local_symbol, MacroTinkerMod, self.scope, self.symbol)
        return line


class ParsedFile:
    def split_by_comma_outside_parentheses(self, s: str) -> List[str]:
        # real*8,dimension(:,:),allocatable -> [real*8, dimension(:,:), allocatable]
        return re.split(r',\s*(?![^()]*\))', s)

    def expressions_within_parentheses(self, s: str) -> str:
        # a -> a; b6(0:a,0:a) -> 0:a,0:a; b5(a,a) -> a,a
        if '(' not in s:
            return s
        else:
            return re.search(r'\((.*?)\)', s).group(1)

    def fortran_double_precision_expression(self, s: str) -> str:
        result = re.search(r'[+-]?(\d+(\.\d*)?|\.\d+)([dDeE][+-]?\d+)?', s)
        if result:
            s2 = result.group(0)
            s2 = s2.replace('d', 'e')
            return s2.replace('D', 'e')
        else:
            return s

    def convert_range_to_length(self, e: str) -> str:
        # : -> :; 100 -> 100; 0:100 -> 101; a:100 -> 101-a; 0:a -> a+1; 100:100+a -> 100+a-99; a:2*a -> 1+2*a-a
        totlen = e
        if ':' == e:
            # :
            pass
        elif ':' in e and len(e) > 1:
            r = e.split(':')

            back, iback, nback = r[1], False, -1
            try:
                nback, iback = int(back), True
            except ValueError:
                pass
            begin, ibegin, nbegin = r[0], False, -1
            try:
                nbegin, ibegin = int(begin), True
            except ValueError:
                pass

            if iback and ibegin:
                # 0:100
                totlen = '{}'.format(nback + 1 - nbegin)
            elif iback and not ibegin:
                # a:100
                totlen = '{}-{}'.format(nback + 1, begin)
            elif not iback and ibegin:
                # 0:a, 100:100+a
                if nbegin == 0:
                    totlen = '{}+1'.format(back)
                else:
                    totlen = '{}-{}'.format(back, nbegin - 1)
            else:
                # a:2*a
                totlen = '1+{}-{}'.format(back, begin)
        else:
            # 100
            pass
        return totlen

    def __init__(self, file_type: FileType, minfo: FileInfo) -> None:
        assert(minfo.file_type == file_type)
        self.current_scope: str = ''
        self.scope_names: List[str] = []
        self.scopes = {}

        if minfo.file_type == FileType.module:
            for line in minfo.raw:
                self.parseModuleLine(line)
        elif minfo.file_type == FileType.subroutine:
            for line in minfo.raw:
                self.parseRoutineLine(line)

    def parseModuleLine(self, raw_line: str) -> None:
        if len(raw_line):
            if raw_line[0] != ' ':
                return

        line: str = raw_line.lstrip()
        if len(line) == 0:
            return

        is_comment: bool = (line[0] == '!' or line[0] == '#')
        word = line.split()[0]
        word1 = word.split(',')[0]
        keys = KnownTypes().keys()
        if is_comment:
            pass
        elif 'implicit ' in line and ' none' in line:
            pass
        elif line == 'save':
            pass
        elif 'module ' == line[0:7]:
            self.current_scope = line.split()[-1]
            m: Variable = Variable()
            m.type = 'module'
            m.symbol = self.current_scope
            self.scope_names.append(self.current_scope)
            self.scopes[self.current_scope] = m
        elif line == 'end':
            non_const: int = 0
            m: Variable = self.scopes[self.current_scope]
            for e in m.entries:
                if not e.is_const:
                    non_const = non_const + 1
            if non_const:
                m.module_has_symbol = True
            self.current_scope = ''
        elif line[0:4] == 'use ':
            m: Variable = self.scopes[self.current_scope]
            m.module_depends_on.append(line.split()[-1])
        elif word in keys or word1 in keys or word in ['parameter']:
            self.parseVariableDefinition(line)
        else:
            print('#error Cannot parse this line: {}'.format(line))

    def parseVariableDefinition(self, line: str) -> None:
        m: Variable = self.scopes[self.current_scope]

        part1: str = ''
        part2: str = ''
        if '::' in line:
            y: List[str] = line.split('::')
            part1, part2 = y[0], y[1]
        else:
            y: List[str] = line.split()
            part1, part2 = y[0], ''.join(y[1:])
        # clear whitespace
        type_aspects: str = ''.join(part1.split(' '))
        symbol_aspects: str = ''.join(part2.split(' '))

        # parse type_aspects
        type_array: List[str] = self.split_by_comma_outside_parentheses(type_aspects)
        this_type: str = ''
        for x in type_array:
            if x in KnownTypes().keys():
                # character*7 -> char
                this_type = KnownTypes()[x]
                break

        this_is_const: bool = False
        if 'parameter' in type_array:
            this_is_const = True

        type_dimension: List[str] = []
        for x in type_array:
            if 'character*' in x:
                # character*7
                y = x.split('*')
                type_dimension.append(y[1])
        for x in type_array:
            if 'dimension(' in x:
                # dimension(2), dimension(:,:)
                ranges: List[str] = self.expressions_within_parentheses(x).split(',')
                for r in ranges:
                    totlen = self.convert_range_to_length(r)
                    type_dimension.append(totlen)

        # parse symbol_aspects
        symbol_array: List[str] = self.split_by_comma_outside_parentheses(symbol_aspects)
        for sa in symbol_array:
            is_an_update: bool = False
            this_dimension: List[str] = [d for d in type_dimension]
            if '=' in sa:
                s: str = self.expressions_within_parentheses(sa)
                tmp: List[str] = s.split('=')
                this_symbol: str = tmp[0]
                this_value: str = tmp[-1]
                this_value = self.fortran_double_precision_expression(this_value)
                assert(this_is_const)
                if this_type != '':
                    # integer, parameter :: x=1000, y=a
                    pass
                else:
                    # parameter (a=1000)
                    is_an_update = True
                    found, e = m.find_entry_by_symbol(this_symbol)
                    if found:
                        e.value = this_value
                        e.is_const = this_is_const
            else:
                # integer a
                # integer b(3,3)
                this_symbol: str = sa.split('(')[0]
                if '(' in sa:
                    ranges: List[str] = self.expressions_within_parentheses(sa).split(',')
                    for r in ranges:
                        totlen: str = self.convert_range_to_length(r)
                        this_dimension.append(totlen)
            if not is_an_update and m.type == 'module':
                e: Variable = Variable()
                e.scope = self.current_scope
                e.symbol = this_symbol
                e.type = this_type
                e.is_const = this_is_const
                e.dimension = this_dimension
                m.entries.append(e)
            elif not is_an_update and m.type == 'fptr':
                found, e = m.find_entry_by_symbol(this_symbol)
                if found:
                    # x: function erfc (x)
                    e.type, e.symbol = this_type, this_symbol
                elif this_symbol == self.current_scope:
                    # erfc: function erfc (x)
                    assert(m.return_type == '')
                    m.return_type = this_type

    def parseRoutineLine(self, raw_line: str) -> None:
        if len(raw_line):
            if raw_line[0] != ' ':
                return

        line: str = raw_line.lstrip()
        if len(line) == 0:
            return

        is_comment: bool = (line[0] == '!' or line[0] == '#')
        word = line.split()[0]
        word1 = word.split(',')[0]
        keys = KnownTypes().keys()
        if is_comment:
            pass
        elif 'implicit ' in line and ' none' in line:
            pass
        elif 'subroutine ' == line[0:11] or 'function ' == line[0:9]:
            e: Variable = Variable()
            vs: List[str] = line.split()
            self.current_scope = vs[1]
            e.symbol = vs[1]
            e.type = 'fptr'
            if vs[0] == 'subroutine':
                e.return_type = 'void'
            args_line: str = self.expressions_within_parentheses(''.join(vs[2:]))
            args_list: List[str] = []
            if args_line != '':
                args_list = args_line.split(',')
            for arg in args_list:
                e1: Variable = Variable()
                e1.symbol = arg
                e.entries.append(e1)
            self.scope_names.append(self.current_scope)
            self.scopes[self.current_scope] = e
        elif word in keys or word1 in keys:
            self.parseVariableDefinition(line)
        elif line == 'end':
            self.current_scope = ''


################################################################################
#                                                                              #
# Writer                                                                       #
#                                                                              #
################################################################################


def printShell(lang: str, files: List[str]) -> None:
    def cmakeContents(files: FileInfoLists) -> str:
        lines: List[str] = []
        relDir = '../source/'
        for f in files.module_files:
            _, stem, ext = getDirnameStemExt(f.filename)
            lines.append(relDir + stem + ext)
        for f in files.subroutine_files:
            _, stem, ext = getDirnameStemExt(f.filename)
            lines.append(relDir + stem + ext)
        contents: str = CMAKE_LIBTINKER.format(FORT_FILES='\n'.join(lines),
            CMAKE_Fortran_COMPILER_ID='{CMAKE_Fortran_COMPILER_ID}', CMAKE_CURRENT_SOURCE_DIR='{CMAKE_CURRENT_SOURCE_DIR}')
        dirname, _, _ = getDirnameStemExt(__file__)
        version_str: str = [line.strip() for line in open(dirname + '/version.txt')][0]
        contents = contents + '\n\n' + '# Generated by v{}.'.format(version_str)
        return contents

    if len(files):
        if lang == 'c':
            print(ShellPart1.format(ModHeaderName='modc'))
        elif lang == 'cpp':
            print(ShellPart1.format(ModHeaderName='modcpp'))
        file_lists: FileInfoLists = FileInfoLists(files)
        for f in file_lists.module_files:
            dirname, stem, _ = getDirnameStemExt(f.filename)
            if lang == 'cpp':
                print('python3 {} --c++ {}/{}.f > $DIR/{}.hh'.format(sys.argv[0], dirname, stem, stem))
            elif lang == 'c':
                print('python3 {} --c {}/{}.f > $DIR/{}.hh'.format(sys.argv[0], dirname, stem, stem))
        for f in file_lists.module_files:
            _, stem, _ = getDirnameStemExt(f.filename)
            print('echo \'#include "detail/{}.hh"\' >> $HEADER'.format(stem))
        if lang == 'cpp':
            print('echo \'#define TINKER_FORTRAN_MODULE_CPP\' > $MODCPP')
            print(R'''echo '#include' '"'$HEADER'"' >> $MODCPP''')
            print('python3 {} --fsrc > routines.cpp'.format(sys.argv[0]))
        elif lang == 'c':
            print('python3 {} --fsrc > routines.c'.format(sys.argv[0]))
        print('python3 {} --func {} > routines.h'.format(sys.argv[0], ' '.join(files)))
        print('echo \'{}\' > ../../CMakeLists.txt'.format(cmakeContents(file_lists)))
        print(ShellPart2)

def printModuleHeaderInCpp(format: str, filenames: List[str]) -> None:
    assert(len(filenames) == 1)
    file_lists: FileInfoLists = FileInfoLists(filenames)
    assert(len(file_lists.module_files))
    mod_list: List[Variable] = []
    for f in file_lists.module_files:
        pmFile: ParsedFile = ParsedFile(FileType.module, f)
        for mk in pmFile.scope_names:
            mod_list.append(pmFile.scopes[mk])
    depends = set()
    for m in mod_list:
        for d in m.module_depends_on:
            depends.add(d)
    depends_list: List[str] = [d for d in depends]
    depends_list.sort()
    lines: List[str] = []
    lines.append('#pragma once')
    lines.append('')
    lines.append('#include "macro.hh"')
    for d in depends_list:
        lines.append('#include "{}.hh"'.format(d))
    print('\n'.join(lines))
    for m in mod_list:
        if format == 'c++':
            print(exportModuleToCpp(m))
        elif format == 'c':
            print(exportModuleToC(m))

def exportModuleToCpp(m: Variable) -> str:
    assert(m.type == 'module')
    lines: List[str] = []
    lines.append('')
    mod_name: str = m.symbol
    if mod_name in CppKeywords():
        mod_name = mod_name + '_'
    lines.append('namespace tinker { namespace %s {' % mod_name)
    if len(m.module_depends_on):
        for d in m.module_depends_on:
            lines.append('using namespace {};'.format(d))
        lines.append('')
    for e in m.entries:
        l = e.f2c(FileType.t9_extern_ref)
        if l != '':
            lines.append(l)
    if m.module_has_symbol:
        lines.append('')
        lines.append('#ifdef TINKER_FORTRAN_MODULE_CPP')
        for e in m.entries:
            l = e.f2c(FileType.t9_extern_c)
            if l != '':
                lines.append(l)
        lines.append('')
        for e in m.entries:
            l = e.f2c(FileType.t9_definition)
            if l != '':
                lines.append(l)
        lines.append('#endif')
    lines.append('} }')
    return '\n'.join(lines)

def exportModuleToC(m: Variable) -> str:
    assert(m.type == 'module')
    lines: List[str] = []
    lines.append('')
    lines.append('''#ifdef __cplusplus
extern "C" {
#endif''')
    for e in m.entries:
        l = e.f2c(FileType.t9_extern_c_with_macro_const)
        if l != '':
            lines.append(l)
    lines.append('''#ifdef __cplusplus
}
#endif''')
    return '\n'.join(lines)

def printFuncHeaderInC(filenames: List[str]) -> None:
    lines: List[str] = []
    file_lists: FileInfoLists = FileInfoLists(filenames)
    length, count = len(file_lists.subroutine_files), 0
    if length:
        lines.append('#pragma once\n')
        lines.append(R'#ifdef __cplusplus')
        lines.append('extern "C" {')
        lines.append(R'#endif')
        lines.append('')
        lines.append(FCharLenTDefinition)
        lines.append(FortranRuntimeHeader)
    for f in file_lists.subroutine_files:
        prFile: ParsedFile = ParsedFile(FileType.subroutine, f)
        if len(prFile.scope_names):
            _, stem, ext = getDirnameStemExt(f.filename)
            lines.append('// {}{}'.format(stem, ext))
        for rk in prFile.scope_names:
            r = prFile.scopes[rk]
            lines.append(exportRoutine(r))
        count = count + 1
        if count != length:
            lines.append('')
    if length:
        lines.append('')
        lines.append(R'#ifdef __cplusplus')
        lines.append('}')
        lines.append(R'#endif')
    print('\n'.join(lines))
    return

def exportRoutine(r: Variable) -> str:
    assert(r.type == 'fptr')
    part1: str = '{} {}_'.format(r.return_type, r.symbol)
    p2_lst: List[str] = []
    char_lst: List[str] = []
    for arg in r.entries:
        if arg.type == 'fptr':
            p2_lst.append(Tinker8.fortranExternalCallInC(r.symbol, arg.symbol))
        else:
            p2_lst.append('{}* {}'.format(arg.type, arg.symbol))
            if arg.type == 'char':
                char_lst.append('tinker_fchar_len_t {}_cap'.format(arg.symbol))
    part2 = ', '.join(p2_lst)
    for ch in char_lst:
        part2 = part2 + ', {}'.format(ch)
    line: str = '{}({});'.format(part1, part2)
    if len(char_lst) == 0:
        line = line + '\n' + '#define tinker_f_{0} {0}_'.format(r.symbol)
    else:
        p3_lst: List[str] = []
        p2_lst.clear()
        char_lst.clear()
        for arg in r.entries:
            if arg.type == 'fptr':
                p3_lst.append(arg.symbol)
                p2_lst.append(Tinker8.fortranExternalCallInC(r.symbol, arg.symbol))
            else:
                if arg.type == 'char':
                    p3_lst.append('{}.string'.format(arg.symbol))
                    char_lst.append(arg.symbol)
                    p2_lst.append('tinker_fchars {}'.format(arg.symbol))
                else:
                    p3_lst.append(arg.symbol)
                    p2_lst.append('{}* {}'.format(arg.type, arg.symbol))
        for ch in char_lst:
            p3_lst.append('{}.capacity'.format(ch))
        call_str = 'return {}_({});'.format(r.symbol, ', '.join(p3_lst))
        line2 = '''inline {} tinker_f_{}({}) {{
    {}
}}'''.format(r.return_type, r.symbol, ', '.join(p2_lst), call_str)
        line = line + '\n' + line2
    return line


################################################################################
#                                                                              #
# main                                                                         #
#                                                                              #
################################################################################


def printHelpMessage() -> None:
    msg: str = '''SYNOPSIS
	{0} [OPTION]... [FILE(S)]...

DESCRIPTION
	Parse the Fortran Tinker source code and generate the C/C++ interface.

	--lang=<LANG> [FILE(S)]...
		print the shell commands for the target language to stdout
		LANG = c cpp
		used as:
			{0} --lang=c [FILE(S)]... | bash
			{0} --lang=cpp [FILE(S)]... | bash

	-h
		print the help message


	--c [FILE]
		print the C header for the Tinker module to stdout
		only accept one file

	--c++ [FILE]
		print the C++ header for the Tinker module to stdout
		only accept one file

	--fsrc
		print the extra source code for the Fortran runtime to stdout

	--func [FILE(S)]...
		print the C header for the Tinker routine(s) to stdout'''.format('parse.py')

    print(msg)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        printHelpMessage()
    else:
        cmd: str = sys.argv[1]
        if cmd == '-h':
            printHelpMessage()
        elif cmd == '--lang=c':
            printShell('c', sys.argv[2:])
        elif cmd == '--lang=cpp':
            printShell('cpp', sys.argv[2:])
        elif cmd == '--c':
            printModuleHeaderInCpp('c', sys.argv[2:])
        elif cmd == '--c++':
            printModuleHeaderInCpp('c++', sys.argv[2:])
        elif cmd == '--fsrc':
            print(FortranRuntimeSrc)
        elif cmd == '--func':
            printFuncHeaderInC(sys.argv[2:])
        else:
            printHelpMessage()
