module X
  class Multivector2D < Multivector
    attr_accessor :e0, :e1, :e2, :e12

    def initialize(input)
      @e0 = input[0]
      @e1 = input[1]
      @e2 = input[2]
      @e12 = input[3]
      @obj_type = "plaintext"
    end

    def to_s
      "#{self.e0}e0 + #{self.e1}e1 + #{self.e2}e2 + #{self.e12}e12"
    end

    def inspect
      to_s
    end

    def geometric_product(m2)
      m = self.clone
      m.e0 = (self.e0*m2.e0) + (self.e1*m2.e1) + (self.e2*m2.e2) - (self.e12*m2.e12)
      m.e1 = (self.e0*m2.e1) + (self.e1*m2.e0) - (self.e2*m2.e12) + (self.e12*m2.e2)
      m.e2 = (self.e0*m2.e2) + (self.e1*m2.e12) + (self.e2*m2.e0) - (self.e12*m2.e1)
      m.e12 = (self.e0*m2.e12) + (self.e1*m2.e2) - (self.e2*m2.e1) + (self.e12*m2.e0)
      m
    end

    def edge_product(m2)
      m = self.clone
      m.e0 = (self.e0*m2.e0) + (self.e1*m2.e1) + (self.e2*m2.e2) + (self.e12*m2.e12)
      m.e1 = (self.e0*m2.e1) + (self.e1*m2.e0) + (self.e2*m2.e12) + (self.e12*m2.e2)
      m.e2 = (self.e0*m2.e2) + (self.e1*m2.e12) + (self.e2*m2.e0) + (self.e12*m2.e1)
      m.e12 = (self.e0*m2.e12) + (self.e1*m2.e2) + (self.e2*m2.e1) + (self.e12*m2.e0)
      m
    end

    alias_method :gp, :geometric_product
    alias_method :ep, :edge_product

    def clifford_conjugation
      signs = [1,-1,-1,-1]
      m = self.clone
      m.data self.data.map.with_index{|d,i| signs[i] == 1 ? d : -d}; m
    end

    alias_method :cc, :clifford_conjugation

    def reverse
      signs = [1,1,1,-1]
      m = self.clone
      m.data self.data.map.with_index{|d,i| signs[i] == 1 ? d : -d}; m
    end

    def negate
      signs = [-1,-1,-1,-1]
      m = self.clone
      m.data self.data.map.with_index{|d,i| signs[i] == 1 ? d : -d}; m
    end

    def amplitude_squared
      self.gp(self.cc)
    end

    def rationalize
      self.amplitude_squared.gp(amplitude_squared.reverse)
    end

    def scalar_mul(scalar)
      m = self.clone
      m.data self.data.map{|d| d * scalar}; m
    end

    def scalar_div(scalar)
      m = self.clone
      m.data self.data.map{|d| d.is_a?(X::Xp) ? d / scalar : Rational(d,scalar)}; m
    end

    def e
      self.class.new [
        Rational(1,2),
        Rational(1,2),
        Rational(1,2),
        -Rational(1,2)
      ]
    end

    def add(m2)
      m = self.clone
      m.data self.data.map.with_index{|d,i| d + m2.data[i]}; m
    end

    def plus(m)
      k = Rational(1,self.number) + Rational(1,m.number)
      self.scalar_mul(k).ep(m)
    end

    def sub(m2)
      m = self.clone
      m.data self.data.map.with_index{|d,i| d - m2.data[i]}; m
    end

    def mul(m2)
      m = self.clone
      m.data self.data.map.with_index{|d,i| d * m2.data[i]}; m
    end

    def div(m2)
      m = self.clone
      m.data self.data.map.with_index{|d,i| d / m2.data[i]}; m
    end

    def numerator
      self.cc.gp(self.amplitude_squared.reverse)
    end

    def denominator
      self.rationalize.e0
    end

    def inverse
      numerator.scalar_div(denominator)
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

    def pown(n)
      self.scalar_mul(self.number ** (n-1))
    end

    def inv
      self.scalar_mul(Rational(1,self.number**2))
    end

    def twop
      self.scalar_mul(Rational(2**self.number,self.number))
    end

    def ntom(n)
      self.scalar_mul(Rational(n**self.number,self.number))
    end

    def rationalize2
      self.gp(self.cc)
    end

    def q
      m = self.clone
      m.e0 = self.e0**2
      m.e1 = -self.e1**2
      m.e2 = -self.e2**2
      m.e12 = -self.e12**2
      m
    end

    def norm
      Math.sqrt(norm_squared)
    end

    def norm_squared
      self.e0**2 + self.e1**2 + self.e2**2 + self.e12**2
    end

  end
end
