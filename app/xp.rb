module X
  class Xp
    attr_accessor :primes, :c, :d, :residue

    def initialize(primes,c,d)
      @primes = primes
      @c = c
      @d = d
      @residue = Residue.new(@primes,@c,@d)
    end

    def to_s
      "{" + residue.residues.join(", ") + "}"
    end

    def inspect
      to_s
    end

    def to_r
      residue.decode
    end

    def to_i
      residue.residues[0]
    end

    def + b
      xp = self.clone
      xp.residue = residue.add(b.residue)
      xp
    end

    def - b
      xp = self.clone
      xp.residue = residue.sub(b.residue)
      xp
    end

    def * b
      xp = self.clone
      xp.residue = residue.mul(b.residue)
      xp
    end

    def / b
      xp = self.clone
      xp.residue = residue.div(b.residue)
      xp
    end

    def -@
      xp = self.clone
      xp.residue = residue.negate
      xp
    end

    def vector
      residue.residues
    end

  end
end
