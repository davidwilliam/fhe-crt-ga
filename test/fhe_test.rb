require "minitest/autorun"
require Dir.pwd + "/x"

class TestFHE < Minitest::Test
  def setup
    @x = X::FHE.new(4,64)
    @p = @x.p
    @q = @x.q

    @n1 = @x.n1
    @n2 = @x.n2
    @n3 = @x.n3
    @n4 = @x.n4

    @x1 = @x.x1
    @x2 = @x.x2
    @x3 = @x.x3
    @x4 = @x.x4

    @k = @x.k

    @bn = @x.bn
  end

  def test_hensel_encoding_decoding
    m1 = 31
    m2 = Rational(17,5)
    m3 = -28
    m4 = Rational(-9,4)

    m1_h = X::FHE.hensel_encoding(m1, @p)
    m2_h = X::FHE.hensel_encoding(m2, @p)
    m3_h = X::FHE.hensel_encoding(m3, @p)
    m4_h = X::FHE.hensel_encoding(m4, @p)

    assert_equal m1, X::FHE.hensel_decoding(m1_h, @p)
    assert_equal m2, X::FHE.hensel_decoding(m2_h, @p)
    assert_equal m3, X::FHE.hensel_decoding(m3_h, @p)
    assert_equal m4, X::FHE.hensel_decoding(m4_h, @p)
  end

  def test_hensel_packing_unpacking
    m1 = 31
    m2 = Rational(17,5)
    m3 = -28
    m4 = Rational(-9,4)

    m1_mat = X::FHE.random_partition(m1,@p)
    m2_mat = X::FHE.random_partition(m2,@p)
    m3_mat = X::FHE.random_partition(m3,@p)
    m4_mat = X::FHE.random_partition(m4,@p)

    mm1 = X::FHE.isomorphic_hensel_packing(m1_mat, @p, @q, @bn)
    mm2 = X::FHE.isomorphic_hensel_packing(m2_mat, @p, @q, @bn)
    mm3 = X::FHE.isomorphic_hensel_packing(m3_mat, @p, @q, @bn)
    mm4 = X::FHE.isomorphic_hensel_packing(m4_mat, @p, @q, @bn)

    assert_equal m1, X::FHE.hensel_unpacking(mm1, @p, @q)
    assert_equal m2, X::FHE.hensel_unpacking(mm2, @p, @q)
    assert_equal m3, X::FHE.hensel_unpacking(mm3, @p, @q)
    assert_equal m4, X::FHE.hensel_unpacking(mm4, @p, @q)
  end

  def test_encryption_decryption
    m1 = X::FHE.random_number(8)
    m2 = X::FHE.random_number(8)

    c1 = @x.encrypt(m1)
    c2 = @x.encrypt(m2)

    c1_d = @x.decrypt(c1)
    c2_d = @x.decrypt(c2)

    assert_equal m1, c1_d
    assert_equal m2, c2_d
  end

  def test_homomorphic_addition
    m1 = X::FHE.random_number(8)
    m2 = X::FHE.random_number(8)

    c1 = @x.encrypt(m1)
    c2 = @x.encrypt(m2)

    res = c1 + c2

    res_d = @x.decrypt(res)

    assert_equal m1 + m2, res_d
  end

  def test_homomorphic_multiplication
    m1 = X::FHE.random_number(8)
    m2 = X::FHE.random_number(8)

    c1 = @x.encrypt(m1)
    c2 = @x.encrypt(m2)

    res = c1 * c2

    res_d = @x.decrypt(res)

    assert_equal m1 * m2, res_d
  end

  def test_ciphertext_rank
    m1 = X::FHE.random_number(8)
    m2 = X::FHE.random_number(8)
    m3 = X::FHE.random_number(8)
    m4 = X::FHE.random_number(8)

    c1 = @x.encrypt(m1)
    c2 = @x.encrypt(m2)
    c3 = @x.encrypt(m3)
    c4 = @x.encrypt(m4)

    c_m_1 = Matrix[
      c1.data,
      c2.data,
      c3.data,
      c4.data
    ]

    assert_equal 4, c_m_1.rank

    m5 = 0
    m6 = 0
    m7 = 0
    m8 = 0

    c5 = @x.encrypt(m5)
    c6 = @x.encrypt(m6)
    c7 = @x.encrypt(m7)
    c8 = @x.encrypt(m8)

    c_m_2 = Matrix[
      c5.data,
      c6.data,
      c7.data,
      c8.data
    ]

    assert_equal 4, c_m_2.rank
  end

  def test_ciphertext_gcd_and_matrix_characteristics
    m1 = 0
    m2 = 0
    m3 = 0
    m4 = 0

    c1 = @x.encrypt(m1)
    c2 = @x.encrypt(m2)
    c3 = @x.encrypt(m3)
    c4 = @x.encrypt(m4)

    c_m = Matrix[
      c1.data,
      c2.data,
      c3.data,
      c4.data
    ]

    assert @x.p != c_m.det
    assert @x.p != c_m.trace

    assert @x.p != c1.e0.gcd(c1.e0)
    assert @x.p != c1.e0.gcd(c1.e2)
    assert @x.p != c1.e0.gcd(c1.e12)
    assert @x.p != c1.e1.gcd(c1.e2)
    assert @x.p != c2.e1.gcd(c1.e12)
    assert @x.p != c3.e2.gcd(c1.e12)
  end

  def test_ciphertext_laplace_expansion
    m1 = 0
    m2 = 0
    m3 = 0
    m4 = 0

    c1 = @x.encrypt(m1)
    c2 = @x.encrypt(m2)
    c3 = @x.encrypt(m3)
    c4 = @x.encrypt(m4)

    c_m = Matrix[
      c1.data,
      c2.data,
      c3.data,
      c4.data
    ]

    assert @x.p != c_m.laplace_expansion(column: 0)
    assert @x.p != c_m.laplace_expansion(column: 1)
    assert @x.p != c_m.laplace_expansion(column: 2)
    assert @x.p != c_m.laplace_expansion(column: 3)
  end

  def test_ciphertext_gcd_difference
    m1 = 3
    m2 = 3
    m3 = 0

    c1 = @x.encrypt(m1)
    c2 = @x.encrypt(m2)
    c3 = @x.encrypt(m3)

    assert @x.p != (c1.number - c2.number).gcd(@x.q)
    assert @x.p != c3.number.gcd(@x.q)
  end

  def test_ciphertext_basic_cpa
    # learning data
    m1 = X::FHE.random_number(8)
    m2 = X::FHE.random_number(8)
    m3 = X::FHE.random_number(8)
    m4 = X::FHE.random_number(8)

    c1 = @x.encrypt(m1)
    c2 = @x.encrypt(m2)
    c3 = @x.encrypt(m3)
    c4 = @x.encrypt(m4)

    c_m_l = Matrix[
      c1.data,
      c2.data,
      c3.data,
      c4.data
    ]

    p_m_l = Matrix[
      [m1],
      [m2],
      [m3],
      [m4]
    ]

    s_m = c_m_l.inverse * p_m_l

    # testing data

    m5 = X::FHE.random_number(8)
    m6 = X::FHE.random_number(8)
    m7 = X::FHE.random_number(8)
    m8 = X::FHE.random_number(8)

    c5 = @x.encrypt(m1)
    c6 = @x.encrypt(m2)
    c7 = @x.encrypt(m3)
    c8 = @x.encrypt(m4)

    c_m_t = Matrix[
      c5.data,
      c6.data,
      c7.data,
      c8.data
    ]

    p_r_t = c_m_t * s_m

    assert m1 != p_r_t[0,0]
    assert m2 != p_r_t[1,0]
    assert m3 != p_r_t[2,0]
    assert m4 != p_r_t[3,0]
  end

  def test_ciphertext_gcd_trace
    m1 = X::FHE.random_number(8)
    m2 = X::FHE.random_number(8)
    m3 = X::FHE.random_number(8)
    m4 = X::FHE.random_number(8)

    m5 = X::FHE.random_number(8)
    m6 = X::FHE.random_number(8)
    m7 = X::FHE.random_number(8)
    m8 = X::FHE.random_number(8)

    c1 = @x.encrypt(m1)
    c2 = @x.encrypt(m2)
    c3 = @x.encrypt(m3)
    c4 = @x.encrypt(m4)

    c5 = @x.encrypt(m5)
    c6 = @x.encrypt(m6)
    c7 = @x.encrypt(m7)
    c8 = @x.encrypt(m8)

    c_m_1 = Matrix[
      c1.data,
      c2.data,
      c3.data,
      c4.data
    ]

    c_m_2 = Matrix[
      c5.data,
      c6.data,
      c7.data,
      c8.data
    ]

    assert @x.p != c_m_1.trace.gcd(c_m_2.trace)
  end

end
