# -*- coding: utf-8 -*-
# Copyright (C) 2014 nineties
# $Id: hchem.py 2014-08-12 10:47:50 nineties $

#= A clone of Tim Hutton's artificial chemistry simulator. =

import numpy as np
import numpy.linalg as la
import pygame, pickle
import re, sys, os, time, math, datetime
import Tkinter as tk
import tkMessageBox

class HChemRule:
    def __init__(self, filename):
        self.cnt = 0
        self.num = None
        self.fill = []
        self.types = []
        self.map = {}
        self.color_count = 0
        self.colors = []
        self.colormap = {}
        self.wildcards = ['X', 'Y']
        self.wildstates = ['x', 'y']
        self.state_max = 10
        self.name = []
        self.ruleb = {} # Rules for bounded pair
        self.ruleu = {} # Rules for unbounded pair
        self.rule_texts = []
        self.colortable = [
            (0, 0, 0), (255, 0, 0), (0, 255, 0), (0, 0, 255),
            (255, 255, 0), (255, 0, 255), (0, 255, 255),
            (255, 127, 0), (255, 0, 127), (127, 255, 0), (127, 0, 255),
            (0, 255, 127), (0, 127, 255)
        ]

        self.parse(filename)
        for t in self.types:
            for s in range(0, self.state_max+1):
                self.to_index(t, str(s))

    def gen_color(self, a):
        self.color_count += 1
        return self.colortable[(self.color_count - 1) % len(self.colortable)]

    def to_index(self, t, n):
        term = t + n
        if not term in self.map:
            self.map[term] = self.cnt
            self.colors.append(self.colormap[t])
            self.name.append(term)
            self.cnt += 1
        return self.map[term]

    def is_valid_type(self, term):
        if term in self.map:
            return True
        else:
            return False

    def get_index(self, term):
        return self.map[term]

    def get_name(self, idx):
        return self.name[idx]

    def parse_expr(self, str):
        if "-" in str:
            bnd = True
            M0, M1 = str.split("-")
        else:
            bnd = False
            M0, M1 = str.split()
        p = re.search(r'([a-wzA-Z]+)(\d+|[xy])', M0)
        q = re.search(r'([a-wzA-Z]+)(\d+|[xy])', M1)
        return (p.group(1), p.group(2), q.group(1), q.group(2), bnd)

    def add_rule(self, L0, l0, L1, l1, lbnd, R0, r0, R1, r1, rbnd, prob):
        if l0 in self.wildstates and l1 in self.wildstates:
            if l0 == l1:
                for s in range(self.state_max+1):
                    s = str(s)
                    _l0 = s
                    _l1 = s
                    if r0 == l0: _r0 = s
                    else:        _r0 = r0
                    if r1 == l0: _r1 = s
                    else:        _r1 = r1
                    self.add_rule(L0, _l0, L1, _l1, lbnd, R0, _r0, R1, _r1, rbnd, prob)
                return
            else:
                for s0 in range(self.state_max+1):
                    for s1 in range(self.state_max+1):
                        if l0 == l1 and s0 > s1: continue
                        s0 = str(s0)
                        s1 = str(s1)
                        _l0 = s0
                        _l1 = s1
                        if r0 == l0:   _r0 = s0
                        elif r0 == l1: _r0 = s1
                        else:          _r0 = r0
                        if r1 == l0:   _r1 = s0
                        elif r1 == l1: _r1 = s1
                        else:          _r1 = r1
                        self.add_rule(L0, _l0, L1, _l1, lbnd, R0, _r0, R1, _r1, rbnd, prob)
                return
        elif l0 in self.wildstates or l1 in self.wildstates:
            if l0 in self.wildstates:
                for s in range(self.state_max + 1):
                    s = str(s)
                    _l0 = s
                    if r0 == l0: _r0 = s
                    else:        _r0 = r0
                    if r1 == l0: _r1 = s
                    else:        _r1 = r1
                    self.add_rule(L0, _l0, L1, l1, lbnd, R0, _r0, R1, _r1, rbnd, prob)
                return
            else:
                for s in range(self.state_max + 1):
                    s = str(s)
                    _l1 = s
                    if r0 == l1: _r0 = s
                    else:        _r0 = r0
                    if r1 == l1: _r1 = s
                    else:        _r1 = r1
                    self.add_rule(L0, l0, L1, _l1, lbnd, R0, _r0, R1, _r1, rbnd, prob)
                return
        LL0 = self.to_index(L0, l0)
        LL1 = self.to_index(L1, l1)
        RR0 = self.to_index(R0, r0)
        RR1 = self.to_index(R1, r1)
        if lbnd:
            if not (LL0, LL1) in self.ruleb:
                self.ruleb[(LL0, LL1)] = []
                if LL0 != LL1:
                    self.ruleb[(LL1, LL0)] = []
            self.ruleb[(LL0, LL1)].append((RR0, RR1, rbnd, prob))
            if LL0 != LL1:
                self.ruleb[(LL1, LL0)].append((RR1, RR0, rbnd, prob))
        else:
            if not (LL0, LL1) in self.ruleu:
                self.ruleu[(LL0, LL1)] = []
                if LL0 != LL1:
                    self.ruleu[(LL1, LL0)] = []
            self.ruleu[(LL0, LL1)].append((RR0, RR1, rbnd, prob))
            if LL0 != LL1:
                self.ruleu[(LL1, LL0)].append((RR1, RR0, rbnd, prob))

    def parse_rule(self, line):
        try:
            lhs, rhs = line.split("->")
        except:
            return
        self.rule_texts.append(line)
        prob = 1.0
        if ":" in rhs:
            rhs, p = rhs.split(":")
            prob = eval(p.strip())
        L0, l0, L1, l1, lbnd = self.parse_expr(lhs.strip())
        R0, r0, R1, r1, rbnd = self.parse_expr(rhs.strip())
        if L0 in self.wildcards and L1 in self.wildcards:
            if L0 == L1:
                for t in self.types:
                    _L0 = t
                    _L1 = t
                    if R0 == L0: _R0 = t
                    else:        _R0 = R0
                    if R1 == L0: _R1 = t
                    else:        _R1 = R1
                    self.add_rule(_L0, l0, _L1, l1, lbnd, _R0, r0, _R1, r1, rbnd, prob)
            else:
                for t0 in self.types:
                    for t1 in self.types:
                        if l0 == l1 and t0 > t1: continue
                        _L0 = t0
                        _L1 = t1
                        if R0 == L0:   _R0 = t0
                        elif R0 == L1: _R0 = t1
                        else:          _R0 = R0
                        if R1 == L0:   _R1 = t0
                        elif R1 == L1: _R1 = t1
                        else:          _R1 = R1
                        self.add_rule(_L0, l0, _L1, l1, lbnd, _R0, r0, _R1, r1, rbnd, prob)
        elif L0 in self.wildcards or L1 in self.wildcards:
            if L0 in self.wildcards:
                for t in self.types:
                    _L0 = t
                    if R0 == L0: _R0 = t
                    else:        _R0 = R0
                    if R1 == L0: _R1 = t
                    else:        _R1 = R1
                    self.add_rule(_L0, l0, L1, l1, lbnd, _R0, r0, _R1, r1, rbnd, prob)
            else:
                for t in self.types:
                    _L1 = t
                    if R0 == L1: _R0 = t
                    else:        _R0 = R0
                    if R1 == L1: _R1 = t
                    else:        _R1 = R1
                    self.add_rule(L0, l0, _L1, l1, lbnd, _R0, r0, _R1, r1, rbnd, prob)
        else:
            self.add_rule(L0, l0, L1, l1, lbnd, R0, r0, R1, r1, rbnd, prob)

    def add_type(self, t):
        self.types.append(t)
        self.colormap[t] = self.gen_color(t)

    def setup_types(self, str):
        lhs, rhs = str.split(":")
        for t in rhs.split(","):
            self.add_type(t.strip())

    def setup_fill(self, str):
        lhs, rhs = str.split(":")
        for decl in rhs.split(","):
            t, p = decl.strip().split(" ")
            self.fill.append((t, eval(p)))

    def parse(self, filename):
        f = open(filename, "r")
        while True:
            line = f.readline()
            if not line: break
            if line[0] == '#': continue
            line = line.strip()
            if line.find("type") == 0:
                self.setup_types(line)
            elif line.find("number of particles") == 0:
                self.num = int(line.split(":")[1].strip())
            elif line.find("state max") == 0:
                self.state_max = int(line.split(":")[1].strip())
            elif line.find("fill") == 0:
                self.setup_fill(line)
            elif "->" in line:
                self.parse_rule(line)
        f.close()

    def check(self, L0, L1, bound):
        if bound:
            if not (L0, L1) in self.ruleb:
                return None
            return self.ruleb[(L0, L1)]
        else:
            if not (L0, L1) in self.ruleu:
                return None
            return self.ruleu[(L0, L1)]

class HChem:
    # n   : number of particles
    # r   : radious of particles
    # v0  : initial velocity of particles
    # dt  : duration of one frame
    # k   : strength of bonds
    # w,h : width and height of the universe
    # seed: random seed
    def __init__(self, filename, n = 1000, r = 10, v0 = None, dt = 0.1,
                 width = 1200, height = 700, bucket_size = None, seed=None):
        self.rule = HChemRule(filename)
        if seed: np.random.seed(seed)
        if v0 == None: v0 = r
        if bucket_size == None: bucket_size = 2*r

        self.n = n
        if self.rule.num: self.n = self.rule.num
        self.r = r
        self.dt = dt
        self.w = width
        self.h = height
        self.speed = 10

        # Initialize positions of particles
        self.pos = np.zeros((n, 2))
        self.pos[:,0] = np.random.uniform(r, width-r, n)
        self.pos[:,1] = np.random.uniform(r, height-r, n)

        # Initialize velocities of particles
        direction = np.random.uniform(0, 2*np.pi, n)
        self.vel = np.zeros((n, 2))
        self.vel[:,0] = v0*np.cos(direction)
        self.vel[:,1] = v0*np.sin(direction) 
        # Initialize types
        self.types = np.zeros(n, dtype=int)
        for k in xrange(self.n):
            p = np.random.uniform(0, 1)
            q = 0
            for (t, r) in self.rule.fill:
                q += r
                if p < q:
                    self.types[k] = self.rule.get_index(t)
                    break
        self.stypes = np.zeros(n, dtype=object)
        for k in xrange(self.n):
            self.stypes[k] = self.rule.get_name(self.types[k])

        # self.bonds[i] == list of indexes of particles which is bound to i.
        self.bonds = np.zeros(n, dtype=object)
        for i in xrange(n): self.bonds[i] = []

        # Initialize buckets for compute_collision detection
        self.bucket_size = bucket_size
        self.nbx = int(math.ceil(float(width)/bucket_size))
        self.nby = int(math.ceil(float(height)/bucket_size))
        self.buckets = np.zeros((self.nbx, self.nby), dtype=object)

    def bucket_index(self, x):
        return (min(max(int(x[0]/self.bucket_size), 0), self.nbx-1),
                min(max(int(x[1]/self.bucket_size), 0), self.nby-1))

    def init_bucket(self):
        for i in xrange(self.nbx):
            for j in xrange(self.nby):
                self.buckets[i,j] = []
        for k in xrange(self.n):
            i, j = self.bucket_index(self.pos[k, :])
            self.buckets[i, j].append(k)

    def add_impulse_from_walls(self):
        r = self.r
        for k in xrange(self.n):
            x = self.pos[k, 0]
            y = self.pos[k, 1]
            vx = self.vel[k, 0]
            vy = self.vel[k, 1]
            if (x < r and vx < 0) or (x > self.w-r and vx > 0):
                self.vel[k, 0] += -2*self.vel[k, 0]
                self.vel[k, 1] += 0
            if (y < r and vy < 0) or (y > self.h-r and vy > 0):
                self.vel[k, 0] += 0
                self.vel[k, 1] += -2*self.vel[k, 1]

    def update_state_of_particle_pair(self, k, l):
        if l in self.bonds[k]:
            # bound pair
            rules = self.rule.check(self.types[k], self.types[l], True)
            if rules:
                for r in rules:
                    p = r[3]
                    if np.random.uniform(0, 1) < p:
                        # print "apply:",
                        # print self.rule.get_name(self.types[k]),
                        # print "-",
                        # print self.rule.get_name(self.types[l]),
                        # print " -> ",
                        # print self.rule.get_name(r[0]),
                        # if r[2]:
                        #     print "-",
                        # else:
                        #     print " ",
                        # print self.rule.get_name(r[1])
                        self.types[k] = r[0]
                        self.types[l] = r[1]
                        if not r[2]:
                            self.bonds[k].remove(l)
                            self.bonds[l].remove(k)
                            return False
                        return True
            return True
        else:
            # unbound pair
            rules = self.rule.check(self.types[k], self.types[l], False)
            if rules:
                for r in rules:
                    p = r[3]
                    if np.random.uniform(0, 1) < p:
                        # print "apply:",
                        # print self.rule.get_name(self.types[k]),
                        # print " ",
                        # print self.rule.get_name(self.types[l]),
                        # print " -> ",
                        # print self.rule.get_name(r[0]),
                        # if r[2]:
                        #     print "-",
                        # else:
                        #     print " ",
                        # print self.rule.get_name(r[1])
                        self.types[k] = r[0]
                        self.types[l] = r[1]
                        if r[2]:
                            self.bonds[k].append(l)
                            self.bonds[l].append(k)
                            return True
                        return False
            return False

    def add_impulse_between_unbound_pair(self, k, l, rx, rv, d2):
        if self.update_state_of_particle_pair(k, l):
            return
        d = math.sqrt(d2)
        n = rx/d
        ldt = -n.dot(rv)
        self.vel[k,:] += ldt*n
        self.vel[l,:] -= ldt*n

    def add_impulse_between_bound_pair(self, k, l, rx, rv, d2):
        d = math.sqrt(d2)
        n = rx/d
        c = rx.dot(rv)
        #ldt = -(2*c + 3*(d2-4*self.r*self.r))/(8*d2)
        #self.vel[k,:] += 2*ldt*rx
        #self.vel[l,:] -= 2*ldt*rx
        if (d < 2*self.r and c < 0) or (d > 2*self.r and c > 0):
            ldt = -n.dot(rv)
            self.vel[k,:] += ldt*n
            self.vel[l,:] -= ldt*n

    def add_impulse_between_particles_sub(self, k, i, j):
        if i < 0 or j < 0 or i >= self.nbx or j >= self.nby: return
        for l in self.buckets[i, j]:
            if k >= l: continue
            rx = self.pos[k,:] - self.pos[l,:]
            rv = self.vel[k,:] - self.vel[l,:]
            if rx.dot(rv) >= 0: continue
            d2 = np.sum(rx*rx)
            if d2 > 4*self.r*self.r: continue
            self.add_impulse_between_unbound_pair(k, l, rx, rv, d2)

    def add_impulse_between_particles(self):
        r = self.r

        # add impulse between unbound pairs
        for k in xrange(self.n):
            i,j = self.bucket_index(self.pos[k, :])
            self.add_impulse_between_particles_sub(k, i-1, j)
            self.add_impulse_between_particles_sub(k, i-1, j-1)
            self.add_impulse_between_particles_sub(k, i-1 ,j+1)
            self.add_impulse_between_particles_sub(k, i, j-1)
            self.add_impulse_between_particles_sub(k, i, j)
            self.add_impulse_between_particles_sub(k, i, j+1)
            self.add_impulse_between_particles_sub(k, i+1, j-1)
            self.add_impulse_between_particles_sub(k, i+1, j)
            self.add_impulse_between_particles_sub(k, i+1, j+1)

    def add_impulse_between_bound_particles(self):
        # add impulse between bound pairs
        for k in xrange(self.n):
            for l in self.bonds[k]:
                if k >= l: continue
                rx = self.pos[k,:] - self.pos[l,:]
                rv = self.vel[k,:] - self.vel[l,:]
                d2 = np.sum(rx*rx)
                self.add_impulse_between_bound_pair(k, l, rx, rv, d2)

    def compute_impulse(self):
        self.init_bucket()
        self.add_impulse_from_walls()
        self.add_impulse_between_particles()
        self.add_impulse_between_bound_particles()
        self.pos += self.vel*self.dt

    def change_speed(self, delta):
        self.speed += delta;
        if self.speed < 0:
            self.speed = 1

    def update(self):
        # Update position
        for k in range(self.speed):
            self.compute_impulse()

    def total_energy(self):
        return np.sum(self.vel*self.vel)/2

    def save(self, fname, type):
        if type == "particles":
            self.save_particles(fname)
        else:
            self.save_rules(fname)

    def save_particles(self, fname):
        try:
            with open(fname, "w") as f:
                f.write(repr(self.n)); f.write("\n")
                f.write(repr(self.dt)); f.write("\n")
                for k in xrange(self.n):
                    f.write(self.stypes[k]); f.write(",")
                    f.write(repr(self.pos[k,0])); f.write(",")
                    f.write(repr(self.pos[k,1])); f.write(",")
                    f.write(repr(self.vel[k,0])); f.write(",")
                    f.write(repr(self.vel[k,1])); f.write("\n")
                for k in xrange(self.n):
                    bonds = " ".join(map(str, self.bonds[k]))
                    f.write(bonds); f.write("\n")
        except Exception as e:
            root = tk.Tk()
            root.withdraw()
            tkMessageBox.showerror('Error', str(e))

    def load_particles(self, fname):
        try:
            with open(fname, "r") as f:
                n     = int(f.readline())
                dt    = float(f.readline())
                pos   = []
                vel   = []
                types = []
                bonds = []
                for k in xrange(n):
                    t, p0, p1, v0, v1 = f.readline().strip().split(",")
                    pos.append((float(p0), float(p1)))
                    vel.append((float(v0), float(v1)))
                    if not self.rule.is_valid_type(t):
                        root = tk.Tk()
                        root.withdraw()
                        tkMessageBox.showerror('Error', 'Unknown type:' + t)
                        return
                    types.append(t)
                for k in xrange(n):
                    bonds.append(map(int, f.readline().strip().split()))
                self.n = n
                self.dt = dt
                self.pos = np.array(pos)
                self.vel = np.array(vel)
                self.stypes = np.array(types, dtype=object)
                self.types = np.array(map(lambda t: self.rule.get_index(t), types), dtype=int)
                self.bonds = bonds

        except Exception as e:
            root = tk.Tk()
            root.withdraw()
            tkMessageBox.showerror('Error', str(e))

    def load(self, fname, type):
        if type == "particles":
            self.load_particles(fname)
        else:
            self.load_rules(fname)

    def start_record(self):
        date = datetime.date.today()
        now  = time.time()
        dirname = str(date) + str(now)
        os.mkdir(dirname)
        self.record_dir = dirname

    def record(self, iteration):
        self.save_particles("%s/iteration-%d.dat" % (self.record_dir, iteration))

class HChemViewer:
    RED   = (255, 0, 0)
    BLUE  = (0, 0, 255)
    WHITE = (255, 255, 255)
    BLACK = (0, 0, 0)
    INFO = [
        "(P) play/pause, (F) stepwise, (R) record, (Q) quit, (T) show/hide particle types, (up) Speed up, (down) Speed down",
        "(drag) move particle, (shift + drag) bind/unbind particles, (double click) change type and state of a particle"
    ]

    def __init__(self, sim, w = None, h = None):
        if w == None: w = sim.w
        if h == None: h = sim.h

        self.sim = sim
        pygame.init()
        self.screen = pygame.display.set_mode((w, h)
            #, pygame.DOUBLEBUF | pygame.FULLSCREEN | pygame.HWSURFACE
            )
        pygame.display.set_caption("Tim Hutton's Artificial Chemistry")
        self.fontsize = 18
        self.font = pygame.font.SysFont(None, self.fontsize)
        info_texts = self.INFO + sim.rule.rule_texts
        self.info = map(lambda text: self.font.render(text, False, self.BLUE),
                info_texts)

        self.speed = 10

        # For events
        self.record        = False
        self.shift         = False
        self.play          = False
        self.stepwise      = False
        self.dragged       = False
        self.which_dragged = None
        self.moving        = False
        self.binding       = False
        self.display_types = False
        self.prev_lclick   = time.time()

    def get_clicked(self):
        for k in xrange(self.sim.n):
            d2 = np.sum((self.sim.pos[k,:] - pygame.mouse.get_pos())**2)
            if d2 < self.sim.r**2:
                return k
                break
        return None

    def ask_particle(self):
        i = self.get_clicked()
        dialog = tk.Tk()
        dialog.title("Enter particle type and state")
        type = tk.StringVar(dialog)
        tk.Label(dialog, text = "type").pack(side=tk.LEFT)
        entry = tk.Entry(dialog, textvariable=type)
        entry.focus_force()
        entry.pack(side=tk.LEFT)
        entry.bind('<Return>', lambda evt: dialog.destroy())
        dialog.mainloop()
        return type.get()

    def ask_file(self, title):
        dialog = tk.Tk()
        dialog.title(title)
        filename = tk.StringVar(dialog)
        savetype = tk.StringVar(dialog)
        savetype.set("particles")
        tk.Radiobutton(dialog, text="Particles", variable=savetype, value="particles").pack(anchor = tk.W)
        tk.Radiobutton(dialog, text="Rules", variable=savetype, value="rules").pack(anchor = tk.W)
        tk.Label(dialog, text = "filename").pack(side=tk.LEFT)
        entry = tk.Entry(dialog, textvariable=filename)
        entry.focus_force()
        entry.pack(side=tk.LEFT)
        entry.bind('<Return>', lambda evt: dialog.destroy())
        dialog.mainloop()
        return filename.get(), savetype.get()

    def check_event(self):
        for event in pygame.event.get():
            if event.type == pygame.KEYDOWN:
                key = pygame.key.get_pressed()
                if key[pygame.K_RSHIFT] or key[pygame.K_LSHIFT]:
                    self.shift = True
                if key[pygame.K_r]:
                    self.record = not self.record
                    if self.record:
                        self.sim.start_record()
                if key[pygame.K_q]:
                    sys.exit()
                if key[pygame.K_p]:
                    self.play = not self.play
                if key[pygame.K_UP]:
                    self.sim.change_speed(1)
                if key[pygame.K_DOWN]:
                    self.sim.change_speed(-1)
                if key[pygame.K_s]:
                    fname, type = self.ask_file("Save configuration")
                    if fname: self.sim.save(fname, type)
                if key[pygame.K_l]:
                    fname, type = self.ask_file("Load configuration")
                    if fname: self.sim.load(fname, type)
                if key[pygame.K_t]:
                    self.display_types = not self.display_types
                if key[pygame.K_f]:
                    self.stepwise = True
                    self.play = True
            if event.type == pygame.KEYUP:
                key = pygame.key.get_pressed()
                if not key[pygame.K_RSHIFT] and not key[pygame.K_LSHIFT]:
                    self.shift = False
            if not self.dragged and event.type == pygame.MOUSEBUTTONDOWN:
                self.play = False
                which_clicked = self.get_clicked()
                clicked = not (which_clicked == None)

                # Detect double click
                t = time.time()
                double_click = False
                if t - self.prev_lclick < 1.0/3:
                    double_click = True
                self.prev_lclick = t

                if clicked and double_click:
                    t = self.ask_particle()
                    try:
                        self.sim.stypes[which_clicked] = t
                        self.sim.types[which_clicked] = self.sim.rule.get_index(t)
                    except:
                        pass
                elif clicked:
                    self.dragged = True
                    self.which_dragged = which_clicked
                    if not self.shift: self.moving = True
                    elif self.shift: self.binding = True
            elif self.dragged and event.type == pygame.MOUSEMOTION:
                if self.moving:
                    self.sim.pos[self.which_dragged,:] = pygame.mouse.get_pos()
            elif self.dragged and event.type == pygame.MOUSEBUTTONUP:
                if self.binding:
                     which_clicked = self.get_clicked()
                     clicked = not (which_clicked == None)
                     if clicked and self.which_dragged != which_clicked and\
                             not (which_clicked in self.sim.bonds[self.which_dragged]):
                         self.sim.bonds[self.which_dragged].append(which_clicked)
                         self.sim.bonds[which_clicked].append(self.which_dragged)
                     elif clicked and self.which_dragged != which_clicked and\
                             (which_clicked in self.sim.bonds[self.which_dragged]):
                         self.sim.bonds[self.which_dragged].remove(which_clicked)
                         self.sim.bonds[which_clicked].remove(self.which_dragged)
                self.moving        = False
                self.binding       = False
                self.dragged       = False
                self.which_dragged = None
            elif event.type == pygame.QUIT:
                sys.exit()

    def loop(self):
        iteration = 0
        screen = self.screen
        while True:
            sim = self.sim
            n   = sim.n
            r   = sim.r
            if self.play:
                iteration += 1
                sim.update()
                if self.record:
                    sim.record(iteration)

            if self.stepwise:
                self.play = False
                self.stepwise = False

            pos = sim.pos

            screen.fill(self.WHITE)
            # Draw particles
            for k in xrange(n):
                pygame.draw.circle(screen, sim.rule.colors[sim.types[k]],
                        (int(pos[k,0]), int(pos[k,1])), r, 1)

            if self.display_types:
                for k in xrange(sim.n):
                    t = sim.rule.get_name(sim.types[k])
                    text = self.font.render(t, False, self.BLACK)
                    rect = text.get_rect()
                    rect.centerx = pos[k,0]
                    rect.centery = pos[k,1]
                    self.screen.blit(text, rect)

            # Draw bonds
            for k in xrange(n):
                for l in sim.bonds[k]:
                    if k >= l: continue
                    pygame.draw.line(screen, self.BLACK, pos[k,:], pos[l,:])

            # Other info
            if self.binding:
                pygame.draw.line(screen, self.BLACK,
                    pos[self.which_dragged,:], pygame.mouse.get_pos())

            y = 10
            for i in self.info:
                self.screen.blit(i, (10, y))
                y += self.fontsize
            text = self.font.render(
                    "time = " + str(iteration*sim.dt),
                    False, self.BLUE)
            self.screen.blit(text, (10, y))
            energy = sim.total_energy()
            text = self.font.render(
                    "energy = " + str(energy),
                    False, self.BLUE)
            self.screen.blit(text, (10, y + self.fontsize))
            self.check_event()
            pygame.display.update()

if __name__ == '__main__':
    if len(sys.argv) >= 2:
        sim = HChem(sys.argv[1])
    else:
        sim = HChem("models/universal_enzyme.txt")
    HChemViewer(sim).loop()
