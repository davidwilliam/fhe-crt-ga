module X
  class Multivector2Dff < Multivector
    attr_accessor :e0, :e1, :e2, :e12, :v1, :v2, :w1, :w2, :v, :w, :modulus

    def self.multivector_for(n)
      while true
        data = Array.new(4){X::FHE.random_number(8)}
        m = new(data)
        if m.rationalize2.e0 == n
          break
        end
      end
      m
    end

    def self.multivector_adds_to(n)
      while true
        data = Array.new(4){X::FHE.random_number(8)}
        m = new(data)
        if m.number == n
          break
        end
      end
      m
    end

    def initialize(input,modulus)
      @e0 = input[0] % modulus
      @e1 = input[1] % modulus
      @e2 = input[2] % modulus
      @e12 = input[3] % modulus
      @obj_type = "plaintext"
      @modulus = modulus
    end

    def to_s
      "#{self.e0}e0 + #{self.e1}e1 + #{self.e2}e2 + #{self.e12}e12"
    end

    def inspect
      to_s
    end

    def mod
      m = self.clone
      m.e0 = self.e0 % modulus
      m.e1 = self.e1 % modulus
      m.e2 = self.e2 % modulus
      m.e12 = self.e12 % modulus
      m
    end

    def geometric_product(m2)
      m = self.clone
      m.e0 = (self.e0*m2.e0) + (self.e1*m2.e1) + (self.e2*m2.e2) - (self.e12*m2.e12)
      m.e1 = (self.e0*m2.e1) + (self.e1*m2.e0) - (self.e2*m2.e12) + (self.e12*m2.e2)
      m.e2 = (self.e0*m2.e2) + (self.e1*m2.e12) + (self.e2*m2.e0) - (self.e12*m2.e1)
      m.e12 = (self.e0*m2.e12) + (self.e1*m2.e2) - (self.e2*m2.e1) + (self.e12*m2.e0)
      m.mod
    end

    def edge_product(m2)
      m = self.clone
      m.e0 = (self.e0*m2.e0) + (self.e1*m2.e1) + (self.e2*m2.e2) + (self.e12*m2.e12)
      m.e1 = (self.e0*m2.e1) + (self.e1*m2.e0) + (self.e2*m2.e12) + (self.e12*m2.e2)
      m.e2 = (self.e0*m2.e2) + (self.e1*m2.e12) + (self.e2*m2.e0) + (self.e12*m2.e1)
      m.e12 = (self.e0*m2.e12) + (self.e1*m2.e2) + (self.e2*m2.e1) + (self.e12*m2.e0)
      m.mod
    end

    def wedge_product(m2)
      m = self.clone
      m.e0 = self.e0
      m.e1 = (self.e0*m2.e1) + (self.e1*m2.e0) + (self.e2*m2.e12) + (self.e12*m2.e2)
      m.e2 = (self.e0*m2.e2) + (self.e1*m2.e12) + (self.e2*m2.e0) + (self.e12*m2.e1)
      m.e12 = (self.e0*m2.e12) + (self.e1*m2.e2) + (self.e2*m2.e1) + (self.e12*m2.e0)
      m.mod
    end

    alias_method :gp, :geometric_product
    alias_method :ep, :edge_product
    alias_method :wp, :wedge_product

    def clifford_conjugation
      signs = [1,-1,-1,-1]
      m = self.clone
      m.data self.data.map.with_index{|d,i| signs[i] == 1 ? d : -d}; m.mod
    end

    alias_method :cc, :clifford_conjugation

    def reverse
      signs = [1,1,1,-1]
      m = self.clone
      m.data self.data.map.with_index{|d,i| signs[i] == 1 ? d : -d}; m.mod
    end

    def negate
      signs = [-1,-1,-1,-1]
      m = self.clone
      m.data self.data.map.with_index{|d,i| signs[i] == 1 ? d : -d}; m.mod
    end

    def amplitude_squared
      self.gp(self.cc)
    end

    def rationalize
      self.amplitude_squared.gp(amplitude_squared.reverse)
    end

    def rationalize2
      self.gp(self.cc)
    end

    def scalar_mul(scalar)
      m = self.clone
      m.data self.data.map{|d| d * scalar}; m.mod
    end

    def scalar_div(scalar)
      m = self.clone
      m.data self.data.map{|d| Rational(d,scalar)}; m
    end

    def add(m2)
      m = self.clone
      m.data self.data.map.with_index{|d,i| d + m2.data[i]}; m.mod
    end

    def sub(m2)
      m = self.clone
      m.data self.data.map.with_index{|d,i| d - m2.data[i]}; m.mod
    end

    def numerator
      self.cc.gp(self.amplitude_squared.reverse)
    end

    def denominator
      self.rationalize.e0
    end

    def inverse
      inv_den = mod_inverse(denominator,modulus)
      numerator.scalar_mul(inv_den)
    end

    def data
      [e0,e1,e2,e12]
    end

    def data(given_data=nil)
      if given_data.nil?
        [e0,e1,e2,e12]
      else
        self.e0 = given_data[0]
        self.e1 = given_data[1]
        self.e2 = given_data[2]
        self.e12 = given_data[3]
        self
      end
    end

    def number
      data.inject(:+)
    end

    def pnumber
      (e0 * e1) + (e2 * e12)
    end

    def dimension
      2
    end

    def bits
      number.to_s(2).size
    end

    def to_matrix
      Matrix[[e0 + e1, e2 + e12],[e2 - e12, e0 - e1]]
    end

    def determinant
      m = to_matrix
      m[0,0]* m[1,1] - m[0,1]*m[1,0]
    end

    def to_multivector(matrix)
      mu = self.clone
      mu.e0 = Rational(1,2) * (matrix[0,0] + matrix[1,1])
      mu.e1 = Rational(1,2) * (matrix[0,0] - matrix[1,1])
      mu.e2 = Rational(1,2) * (matrix[0,1] + matrix[1,0])
      mu.e12 = Rational(1,2) * (matrix[0,1] - matrix[1,0])
      mu
    end

    def inverse_m
      m = to_matrix
      m_inv = Rational(1,determinant) * Matrix[[m[1,1],-m[0,1]],[-m[1,0],m[0,0]]]
      to_multivector(m_inv)
    end

    def rank
      to_matrix.rank
    end

    def blade_0
      m = self.clone
      m.e0 = self.e0
      m.e1 = 0
      m.e2 = 0
      m.e12 = 0
      m
    end

    def dual
      scalar_mul(Rational(1,e12))
    end

    def + m2
      self.add(m2)
    end

    def - m2
      self.sub(m2)
    end

    def * m2
      self.gp(m2).scalar_mul(4)
    end

    def lambda1
      e0 + Math.sqrt(e1**2 - e12**2 + e2**2)
    end

    def lambda2
      e0 - Math.sqrt(e1**2 - e12**2 + e2**2)
    end

    def eigenvalues
      [lambda1, lambda2]
    end

    def a
      # supporting multivector for eigenvector1 (v1)
      self.class.new [
        Rational(1,2) * ((e0 + e1 - lambda1) + (e0 - e1 - lambda1)),
        Rational(1,2) * ((e0 + e1 - lambda1) - (e0 - e1 - lambda1)),
        Rational(1,2) * ((e2 + e12) + (e2 + e12)),
        Rational(1,2) * ((e2 + e12) - (e2 + e12))
      ]
    end

    def b
      # supporting multivector for eigenvector2 (w2)
      self.class.new [
        Rational(1,2) * ((e0 + e1 - lambda2) + (e0 - e1 - lambda2)),
        Rational(1,2) * ((e0 + e1 - lambda2) - (e0 - e1 - lambda2)),
        Rational(1,2) * ((e2 + e12) + (e2 + e12)),
        Rational(1,2) * ((e2 + e12) - (e2 + e12))
      ]
    end

    def get_eigenvectors
      v1_ = (a.e2 + a.e12).to_r
      v2_ = Rational(-(a.e0 + a.e1) * v1_, a.e2 + a.e12)
      w1_ = (b.e2 + b.e12).to_r
      w2_ = -Rational((b.e0 + b.e1) * w1_, b.e2 + b.e12)

      # puts "v1_ = #{v1_}, v2_ = #{v2_}"
      # puts "w1_ = #{w1_}, w2_ = #{w2_}"

      v1_gcd_v2 = (v1_.numerator * v2_.denominator).gcd(v2_.numerator * v1_.denominator)

      # puts "v1_gcd_v2 = #{v1_gcd_v2}"

      if v1_gcd_v2 != 1
        @v1 = Rational(v1_,v1_gcd_v2)
        @v2 = Rational(v2_,v1_gcd_v2)
      else
        @v1 = v1_
        @v2 = v2_
      end

      w1_gcd_w2 = (w1_.numerator * w2_.denominator).gcd(w2_.numerator * w1_.denominator)

      # puts "w1_gcd_w2 = #{w1_gcd_w2}"

      if w1_gcd_w2 != 1
        @w1 = Rational(w1_,w1_gcd_w2)
        @w2 = Rational(w2_,w1_gcd_w2)
      else
        @w1 = w1_
        @w2 = w2_
      end

      @v = [v1,v2]
      @w = [w1,w2]

      [@v,@w]
    end

    def p
      get_eigenvectors
      # supporting multivector for eigenvector2 (w2)
      # puts "v1 = #{v1}, v2 = #{v2}"
      # puts "w1 = #{w1}, w2 = #{w2}"
      # puts "eigenvectors = #{eigenvectors}"
      self.class.new [
        Rational(1,2) * (v1 + w2),
        Rational(1,2) * (v1 - w2),
        Rational(1,2) * (w1 + v2),
        Rational(1,2) * (w1 - v2)
      ]
    end

    def d
      p.inverse.gp(self).gp(p)
    end

    def d_pow(n)
      self.class.new [
        Rational(1,2) * ((d.e0 + d.e1)**2 + (d.e0 - d.e1)**2),
        Rational(1,2) * ((d.e0 + d.e1)**2 - (d.e0 - d.e1)**2),
        Rational(1,2) * ((d.e2 + d.e12) + (d.e2 - d.e12)),
        Rational(1,2) * ((d.e2 + d.e12) - (d.e2 - d.e12))
      ]
    end

    def pow(n)
      p.gp(d_pow(n)).gp(p.inverse)
    end

    def e
      self.class.new [
        Rational(1,2),
        Rational(1,2),
        Rational(1,2),
        -Rational(1,2)
      ]
    end

    def pow2(n)
      e_ = self.class.new [
        Rational(1,2) * ((Rational(n,2) + Rational(n,2)) + (Rational(n,2) - Rational(n,2))),
        Rational(1,2) * ((Rational(n,2) + Rational(n,2)) - (Rational(n,2) - Rational(n,2))),
        Rational(1,2) * ((Rational(n,2) + Rational(n,2)) + (Rational(n,2) - Rational(n,2))),
        Rational(1,2) * ((Rational(n,2) + Rational(n,2)) - (Rational(n,2) - Rational(n,2)))
      ]

      puts "e_ = #{e_}"
      self.gp(e_)
    end

    # def mod_inverse(a, m)
    #   a = a % m
    #   for x in (1..m) do
    #     if ((a * x) % m == 1)
    #       return x
    #     end
    #   end
    #   return 1
    # end

    def mod_inverse(num, mod)
      g, a, b = Tools.extended_gcd(num, mod)
      unless g == 1
        raise ZeroDivisionError.new("#{num} has no inverse modulo #{mod}")
      end
      a % mod
    end

    def to_n(n)
      self.scalar_mul(self.number ** (n - 1)).mod
    end

    def to_m(m)
      if m.rationalize.e0.to_i.even?
        self.scalar_mul(self.rationalize2.e0 ** ((m.rationalize2.e0 - 1)/2).to_i).gp(self)
      else
        self.scalar_mul(self.rationalize2.e0 ** ((m.rationalize2.e0 - 1)/2).to_i)
      end
    end

    def bnumber
      data.join.to_i(2)
    end

    def tnumber
      e0 * e1 + e2 * e12
    end

  end
end
