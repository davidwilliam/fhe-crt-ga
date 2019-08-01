module X
  class Residue
    ############################## CLASS MEMBERS ##############################
    attr_accessor  :primes, :residues

    ############################## CONSTRUCTOR ##############################
    # c = numerator
    # d = denominator
    def initialize(primes,c=0,d=1)
      @primes = primes
      @residues = []
      if(c!=0)
        primes.each do |p|
          @residues << Padic.hensel_code(p,1,c,d)
        end
      else
        @residues = Array.new(primes.size,0)
      end
    end

    ############################## OBJECT METHODS ##############################
    def decode_helper
      p = 1
      (0..primes.size-1).each do |i|
        p = p * primes[i]
      end

      # Compute residue
      i = 0
      residue = 0
      residues.each do |r|
        p_factor = p / primes[i]
        p_prime = Tools.mod_inverse(p_factor,primes[i])
        residue = residue + ((p_factor * p_prime * r))
        i+=1
      end

      residue = residue % p

      [p,residue]
    end


    def decode
      p,residue = decode_helper
      Tools.ctr(p, residue , Math.sqrt(p))
    end


    def add(b)
      result = []
      (0..self.primes.size-1).each do |i|
        result << self.residues[i] + b.residues[i]
      end

      c = Residue.new(self.primes)
      c.residues = result
      c
    end

    def sub(b)
      result = []
      (0..self.primes.size-1).each do |i|
        result << self.residues[i] - b.residues[i]
      end

      c = Residue.new(self.primes)
      c.residues = result
      c
    end

    def mul(b)
      result = []
      (0..self.primes.size-1).each do |i|
        result << self.residues[i] * b.residues[i]
      end

      c = Residue.new(self.primes)
      c.residues = result
      c
    end

    def div(b)
      self.mul(b.inverse)
    end

    def inverse
      p,residue = decode_helper
      inv_hensel_code = Tools.mod_inverse(residue,p)
      c = Residue.new(primes)
      c.residues = []
      primes.each do |p|
        c.residues << inv_hensel_code % p
      end
      c
    end



    def negate
      self.mul(Residue.new(self.primes,-1))
    end

    # Prints value
    def to_s
      residues.to_s
    end

    def inspect
      to_s
    end


  end

end
