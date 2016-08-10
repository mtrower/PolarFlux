from uncertainty import Measurement as M

def test_addition_commutativity():
    x = M(3.04, .23)
    y = M(2.96, .14)
    l = x + y
    r = y + x
    assert l.v == r.v and l.u == r.u

def test_multiplication_commutativity():
    x = M(3.04, .23)
    y = M(2.96, .14)
    l = x*y
    r = y*x
    assert l.v == r.v and l.u == r.u

def test_multiple_ops():
    x = M(3.04, .23)
    y = M(2.96, .14)
    z = M(3.00, .02)
    result1 = x*y*z
    result2 = y*z*x
    result3 = z*y*x
    print(result1)
    print(result2)
    print(result3)
    assert result1.v==result2.v and result2.v==result3.v and result1.v==result3.v

def test_array_addition():
    