import re
from copy import deepcopy

file_path = './tests/4.txt'

lst_of_variables = dict()

class Failure(Exception):
    pass

class TempMultEq():
    def __init__(self, S = None, M = None):
        self.S = S
        self.M = M

class Variable():
    def __init__(self, name = None, M = None):
        self.name = name
        self.M = M

class MultiTerm():
    def __init__(self, fsymb = None, args = None):
        self.fsymb = fsymb
        self.args = args
    
class MultiEquation():
    def __init__(self, counter = None, S = None, M = None, var_number = None):
        #текущий счётчик
        self.counter = counter
        #S-часть мультиуравнения
        self.S = S
        #M-часть мультиуравнения
        self.M = M
        #количество переменных в S-части *используется на этапе компактификации*
        self.var_number = var_number

class U():
    def __init__(self, equations = None, zero_counter_multeq = None, multeq_number = None):
        #мультиуравнения
        self.equations = equations
        #мультиуравнения с нулевым счётчиком
        self.zero_counter_multeq = zero_counter_multeq
        #количество мультиуравнений в U-части системы
        self.multeq_number = multeq_number

class System():
    def __init__(self, U = None, T = None):
        #U-Part of a system
        self.U = U
        #T-Part of a system
        self.T = T

temp = r'\s*constructors\s*=\s*(?P<constructors>(?:([A-Za-z]\(\d+\),\s*)*([A-Za-z]\(\d\))))\s*variables\s*=\s*(?P<variables>(?:([A-Za-z],\s*)*[A-Za-z]))\s*first_term\s*=\s*(?P<first_term>.*)\s*second_term\s*=\s*(?P<second_term>.*)'

def parse_data(file_path):

    f = open(file_path, 'r')
    data = f.read().replace('\n', ' ')
    f.close()

    data = re.sub(' +', ' ', data)
    #print(data)

    try:
        match = re.match(temp, data)
    except:
        raise Failure('ошибка при парсинге входных данных')

    constructors = dict()
    constructors.update([(elem[0], int(re.findall('\d+', elem)[0])) for elem in match.group('constructors').replace(' ', '').split(',')])
    #print(f'constructors: {constructors}')

    variables = dict()    
    variables.update([(elem[0], 0) for elem in match.group('variables').replace(' ', '').split(',')])
    #print(f'variables:{variables}')

    first_term = match.group('first_term')
    #print(f'first_term:{first_term}')

    second_term = match.group('second_term')
    #print(f'second_term:{second_term}')

    return constructors, variables, first_term, second_term

def merge_multiterms(first_M, second_M):
    if first_M is None:
        first_M = second_M
    else:
        if not(second_M is None):
            if not(first_M.fsymb == second_M.fsymb):
                raise Failure('не унифицируется')
            else:
                first_arg = first_M.args
                second_arg = second_M.args
                for i in range (len(first_arg)):
                    first_arg[i].S += second_arg[i].S
                    first_arg[i].M = merge_multiterms(first_arg[i].M, second_arg[i].M)
    return first_M

def reducee(M, frontier):
    arg = M.args
    for i in range (len(arg)):
        if arg[i].S == []:
            arg[i].M, frontier = reducee(arg[i].M, frontier)
        else:
            frontier.append(deepcopy(arg[i]))
            arg[i] = TempMultEq(arg[i].S[0], None)
    return M, frontier

def unify(R):
    while R.U.multeq_number != 0:
        if R.U.zero_counter_multeq == []:
            raise Failure('не унифицируется')
        multeq = R.U.zero_counter_multeq[0]
        R.U.zero_counter_multeq.pop(0)
        R.U.multeq_number -= 1
        if multeq.M is not None:
            frontier = []
            multeq.M, frontier = reducee(multeq.M, frontier)
            while frontier != []:
                varss = frontier[0].S
                V = varss[0]
                varss.pop(0)
                V = lst_of_variables[V]
                mult = V.M
                mult.counter -= 1
                while varss != []:
                    V = varss[0]
                    varss.pop(0)
                    V = lst_of_variables[V]
                    first_mult = V.M
                    first_mult.counter -= 1
                    #merge_multeqs
                    if not(mult == first_mult):
                        if mult.var_number < first_mult.var_number:
                            multt = mult
                            mult = first_mult
                            first_mult = multt
                        mult.counter += first_mult.counter
                        mult.var_number += first_mult.var_number
                        varsss = first_mult.S
                        while varsss != []:
                            V1 = varsss[0]
                            varsss.pop(0)
                            V1.M = mult
                            mult.S.append(V1)
                        mult.M = merge_multiterms(mult.M, first_mult.M)
                        R.U.multeq_number -= 1
                        #end_merge_multeqs
                mult.M = merge_multiterms(mult.M, frontier[0].M)
                if mult.counter == 0:
                    R.U.zero_counter_multeq.insert(0, mult)
                frontier.pop(0)
        R.T.append(multeq)
    return R

def find_arg(str, index, variables, constructors):
    for i, ch in enumerate(str):
        if ch in variables or ch in constructors:
            index += i
            return ch, index

def skip_constructor_args(str, index):
    start_index = index
    f = 1
    while f > 0:
        if str[index] == '(':
            f += 1
        elif str[index] == ')':
            f -= 1
        index += 1
    return str[start_index - 1:index], index

def parse_term(term, constructors, variables):
    name = term[0]
    if name in constructors:
        arity = constructors[name]
        return MultiTerm(name, parse_subterms(term[1:], arity, constructors, variables))
    elif name in variables:
        return MultiTerm(None, parse_subterms(term , 1, constructors, variables))
    else:
        raise Failure('ошибка при парсинге терма')

def parse_subterms(subterms, arity, constructors, variables):
    global lst_of_variables
    args = []
    index = 0
    for _ in range(arity):
        arg, index = find_arg(subterms[index:], index, variables, constructors)
        #print(arg, index)
        if arg in variables:
            index += 1
            variables[arg] += 1
            args.append(TempMultEq([arg], None))
        elif arg in constructors:
            arity = constructors[arg]
            if not(arity):
                args.append(TempMultEq([], MultiTerm(arg, [])))
                index += 1
            else:
                new_subterms, new_index = skip_constructor_args(subterms[index:], 2)
                new_index += index
                #print(index, new_index, new_subterms)
                args.append(TempMultEq([], MultiTerm(arg, parse_subterms(new_subterms, arity, constructors, variables))))
                index = new_index
    return args

def form_args(M):
    result = M.fsymb + '('
    for arg in M.args:
        if arg.S != []: result += str(arg.S) + ', '
        else:  
            if arg.M.args != []: result += form_args(arg.M) + ', '
            else: result += str(arg.M.fsymb) + ', '
    result = result[:-2] if M.args != [] else result
    result += ')'
    return result

def main():
    global lst_of_variables
    
    constructors, variables, first_term, second_term = parse_data(file_path)

    for variable in variables:
        lst_of_variables[variable] = Variable(variable)

    first_M = parse_term(first_term, constructors, variables)
    second_M = parse_term(second_term, constructors, variables)
    if not(first_M.fsymb) and second_M.fsymb and not(len(second_M.args)):
        print(f'{first_M.args[0].S[0]} ::= {{{second_M.fsymb}}}')
        return
    elif first_M.fsymb and not(second_M.fsymb):
        print(f'{second_M.args[0].S[0]} ::= {{{first_M.fsymb}}}')
        return  

    M = merge_multiterms(first_M, second_M)

    multeq_number = len(variables) + 1
    equations = [MultiEquation(0, [Variable('x0')], M, 1)]
    equations[0].S[0].M = equations[0]
    for variable in variables:
        if variables[variable]:
            equation = [MultiEquation(variables[variable], [lst_of_variables[variable]], None, 1)]
            equation[0].S[0].M = equation[0]
            equations += equation
        else:
            multeq_number -= 1
    zero_counter_multeq = [equations[0]]
    R = System(U(equations, zero_counter_multeq, multeq_number), [])
    R = unify(R)
    for multeq in R.T:
        vars = [multeq.S[i].name for i in range(len(multeq.S))]
        if multeq.M is None:
            args = [var.name for var in multeq.S]
            print(f'{{{", ".join(args)}}} ::= {{}}')
            continue
        args = form_args(multeq.M)
        print(f'{{{", ".join(vars)}}} ::= {{{args}}}') if not(constructors.get(multeq.M.fsymb))\
        else print(f'{{{", ".join(vars)}}} ::= {{{args}}}')

if __name__ == '__main__':
    main()