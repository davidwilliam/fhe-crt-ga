module X
  class FHE
    attr_accessor :p, :n1, :n2, :n3, :n4,
                  :bn1, :bn2, :bn3, :bn4,
                  :x1, :x2, :x3, :x4,
                  :q, :k, :e, :bn

    def self.random_number(bits)
      OpenSSL::BN::rand(bits).to_i
    end

    def self.random_prime(bits)
      OpenSSL::BN::generate_prime(bits).to_i
    end

    def initialize(gamma,lambda)
      p_bits = gamma * lambda
      q_bits = 2 * gamma * lambda
      ns_bits = 16 * gamma * lambda
      k1_bits = lambda / 2
      k2_bits = lambda / 2

      @p = FHE.random_prime(p_bits)
      @q = FHE.random_prime(q_bits) * FHE.random_prime(q_bits)

      ns = [ FHE.random_number(ns_bits) ]

      while ns.size < 4
        nn = FHE.random_number(ns_bits)
        if ns.map{|n| n.gcd(nn)}.uniq == [1]
          ns << nn
        end
      end

      @n1 = ns[0]; @n2 = ns[1]; @n3 = ns[2]; @n4 = ns[3]

      @bn = @p * @n1 * @n2 * @n3 * @n4

      @bn1 = n2 * n3 * n4
      @bn2 = n1 * n3 * n4
      @bn3 = n1 * n2 * n4
      @bn4 = n1 * n2 * n3

      @x1 = self.class.mod_inverse(bn1, n1)
      @x2 = self.class.mod_inverse(bn2, n2)
      @x3 = self.class.mod_inverse(bn3, n3)
      @x4 = self.class.mod_inverse(bn4, n4)

      valid_key = false
      while valid_key == false
        begin
          k1_ = self.class.random_number(k1_bits)
          k2_ = self.class.random_number(k2_bits)

          k1_mat = self.class.random_partition(k1_,@p)
          k2_mat = self.class.random_partition(k2_,@p)

          k1 = self.class.hensel_packing(k1_mat,@p,@bn)
          k2 = self.class.hensel_packing(k2_mat,@p,@bn)

          @k = k1.gp(k2)

          @k.inverse

          valid_key = true
        rescue => e
          # puts e
        end
      end
    end

    def self.hensel_encoding(n, p)
      Xp.new([p], n.numerator, n.denominator).to_i
    end

    def self.hensel_decoding(n, p)
      Xp.new([p], n.numerator, n.denominator).to_r
    end

    def self.random_partition(m,p)
      p_bits = p.bit_length

      rs = Array.new(4){ FHE.random_number(p_bits / 8) }

      d1 = Rational(rs[0],rs[1])
      d2 = Rational(rs[2],rs[3])

      data = Array.new(3){
        [0,1].sample == 1 ? 1 * FHE.random_number(p_bits / 8) : (-1) * FHE.random_number(p_bits / 8)
      }

      data << (m - data.inject(:+))

      # matrix m
      m_mat = [
        data[0] + d1,
        data[1] - d1,
        data[2] + d2,
        data[3] - d2
      ]

      m_mat
    end

    def self.isomorphic_multivector_packing(p2,q)
      data = [
        Rational(1,2),
        Rational(1,2),
        Rational(1,2),
        Rational(-1,2)
      ]

      data_p = data.map{|d|
        dxp = Xp.new([p2], d.numerator, d.denominator); dxp.to_i
      }

      Multivector2Dff.new data_p, q
    end

    def self.hensel_packing(m_mat,p,bn)
      p_bits = p.bit_length

      data_p = m_mat.map{|d|
        hensel_encoding(d, p)
      }

      mm_p = Multivector2Dff.new data_p, bn

      mm_p
    end

    def self.isomorphic_hensel_packing(m_mat,p,q,bn)
      p_bits = p.bit_length

      a = Multivector2D.new m_mat
      e = Multivector2D.new [
        Rational(1,2),
        Rational(1,2),
        Rational(1,2),
        Rational(-1,2)
      ]
      b = a.gp(e)

      data_p = b.data.map{|d|
        random_number(p_bits / 4) * q + hensel_encoding(d, p)
      }

      mm_p = Multivector2Dff.new data_p, bn

      mm_p
    end

    def self.hensel_unpacking(mm, p, q)
      data_ = mm.data.map{|d|
        hensel_decoding(d % q, p)
      }

      data_.inject(:+)
    end

    def crt_encoding(mm)
      mm_crt = Multivector2Dff.new [
        (mm.e0 * bn1 * x1) % bn,
        (mm.e1 * bn2 * x2) % bn,
        (mm.e2 * bn3 * x3) % bn,
        (mm.e12 * bn4 * x4) % bn
      ], bn

      mm_crt
    end

    def crt_decoding(mm_crt)
      mm = Multivector2Dff.new [
        mm_crt.e0 % n1,
        mm_crt.e1 % n2,
        mm_crt.e2 % n3,
        mm_crt.e12 % n4
      ], bn

      mm
    end

    def self.modular_pow( base, power, mod )
      res = 1
      while power > 0
        res = (res * base) % mod if power & 1 == 1
        base = base ** 2 % mod
        power >>= 1
      end
      res
    end

    def self.mod_inverse(num, mod)
      g, a, b = extended_gcd(num, mod)
      unless g == 1
        raise ZeroDivisionError.new("#{num} has no inverse modulo #{mod}")
      end
      a % mod
    end

    def self.extended_gcd(x, y)
      if x < 0
        g, a, b = extended_gcd(-x, y)
        return [g, -a, b]
      end
      if y < 0
        g, a, b = extended_gcd(x, -y)
        return [g, a, -b]
      end
      r0, r1 = x, y
      a0 = b1 = 1
      a1 = b0 = 0
      until r1.zero?
        q = r0 / r1
        r0, r1 = r1, r0 - q*r1
        a0, a1 = a1, a0 - q*a1
        b0, b1 = b1, b0 - q*b1
      end
      [r0, a0, b0]
    end

    def gp_encrypt(mm)
      k.gp(mm).gp(k.inverse)
    end

    def gp_decrypt(c)
      k.inverse.gp(c).gp(k)
    end

    def encrypt(m)
      m_mat = self.class.random_partition(m,p)
      mm = self.class.isomorphic_hensel_packing(m_mat,p,q,bn)
      c_crt = crt_encoding(mm)
      c = gp_encrypt(c_crt)

      c
    end

    def decrypt(c)
      c_crt = gp_decrypt(c)
      mm = crt_decoding(c_crt)
      m = self.class.hensel_unpacking(mm, p, q)
      m
    end

    def add(c1,c2)
      c1.add(c2)
    end

    def mul(c1,c2)
      c1.gp(c2).scalar_mul(4)
    end

  end
end
