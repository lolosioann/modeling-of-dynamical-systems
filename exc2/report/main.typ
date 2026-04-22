#import "@preview/cetz:0.3.4": canvas
#import "@preview/cetz:0.3.4": draw as cdraw

#set text(lang: "el", size: 11pt)
#set page(
  margin: (x: 2.5cm, y: 2.5cm),
  numbering: "1",
  number-align: center,
  header: context {
    if counter(page).get().first() > 1 [
      #set text(size: 9pt, fill: luma(120))
      #grid(
        columns: (1fr, 1fr),
        align(left)[Προσομοίωση και Μοντελοποίηση — Εργασία 2], align(right)[Λώλος Ιωάννης · 10674],
      )
      #line(length: 100%, stroke: 0.4pt + luma(180))
    ]
  },
)
#set par(justify: true)
#set math.equation(numbering: "(1)")

#align(center)[
  #v(2em)
  #text(size: 16pt)[
    Προσομοίωση και Μοντελοποίηση Δυναμικών Συστημάτων \
    Εργασία 2
  ]
  #v(1.2em)
  #block(width: 65%)[
    #grid(
      columns: (1fr, 1fr),
      align(left)[#text(size: 12pt)[Λώλος Ιωάννης - 10674]], align(right)[#text(size: 10pt, font: "Courier Prime")[lolosioann\@ece.auth.gr]],
    )
  ]
  #v(0.4em)
  #text(size: 11pt)[23 Απριλίου 2026]
  #v(2em)
]

= Εισαγωγή

Η παρούσα αναφορά αφορά την Εκτίμηση Άγνωστων
Παραμέτρων σε πραγματικό χρόνο εφαρμοσμένη σε συστήματα με, ή χωρίς εξωτερικές διαταραχές.

Στο *Θέμα 1* εκτιμούνται οι παράμετροι που περιγράφουν την περιστροφική κίνηση ενός εναέριου οχήματος, με την Μέθοδο Κλίσης.

Στο *Θέμα 2* εφαρμόζεται η Μέθοδος Lyapunov για την on-line εκτίμηση παραμέτρων ενός απλού εκκρεμούς, παρουσία θορύβου.

= Θέμα 1

== Περιγραφή του Συστήματος

Θεωρούμε τη γραμμικοποιημένη εξίσωση προσανατολισμού ενός εναέριου οχήματος γύρω από έναν άξονα περιστροφής:

$ J dot.double(phi)(t) = -k dot(phi)(t) + b u(t) + d(t) $ <eq-plant>

όπου $phi(t)$ [rad] είναι η γωνία προσανατολισμού, $dot(phi)(t)$, $dot.double(phi)(t)$ η γωνιακή ταχύτητα και επιτάχυνση αντίστοιχα, και $u(t)$ [N·m] η ροπή ελέγχου. Η ροπή αδράνειας $J = 0.025$ kg·m² είναι γνωστή, ενώ ο συντελεστής απόσβεσης $k > 0$ και η σταθερά $b > 0$ είναι άγνωστες παράμετροι που επιθυμούμε να εκτιμήσουμε σε πραγματικό χρόνο. Το $d(t)$ αντιπροσωπεύει άγνωστη εξωτερική διαταραχή. Για τα πειράματα επιλέγεται $u(t) = 0.25 sin(0.5 pi t)$ και υποθέτουμε ότι τα σήματα $phi(t)$, $dot(phi)(t)$, $dot.double(phi)(t)$, $u(t)$ είναι μετρήσιμα.

#figure(
  canvas(length: 1cm, {
    import cdraw: *

    // Rotation axis (dashed horizontal line through origin)
    line((-4.2, 0), (4.2, 0), stroke: (dash: "dashed", thickness: 0.5pt, paint: luma(120)))
    content((4.55, 0), text(size: 8pt, fill: luma(90))[άξονας])

    // Tilted vehicle body — rotated by angle phi ~ 20deg
    // Central hub
    circle((0, 0), radius: 0.18, fill: black)

    // Arms of the quadrotor (tilted)
    // Use rotation by ~20 degrees: cos=0.94, sin=0.34
    let c = 0.94
    let s = 0.34
    let L = 2.2
    // Main arm endpoints
    let p1 = (-L * c, -L * s)
    let p2 = (L * c, L * s)
    // Perpendicular arm endpoints (shorter, for visual)
    let q1 = (L * s * 0.55, -L * c * 0.55)
    let q2 = (-L * s * 0.55, L * c * 0.55)

    line(p1, p2, stroke: 2pt)
    line(q1, q2, stroke: 1.2pt, paint: luma(140))

    // Rotors (circles at arm tips)
    circle(p1, radius: 0.32, stroke: 1.2pt, fill: luma(240))
    circle(p2, radius: 0.32, stroke: 1.2pt, fill: luma(240))

    // Reference horizontal line (from hub to the right)
    // line((0, 0), (1.9, 0), stroke: 0.8pt)

    // Angle arc phi between horizontal and main arm
    arc((1.2, 0), start: 0deg, stop: 20deg, radius: 1.3, stroke: 0.9pt)
    content((1.6, 0.28), text(size: 10pt)[$phi(t)$])

    // Torque u(t) — curved arrow around hub
    arc((-0.5, 0.9), start: 120deg, stop: 240deg, radius: 0.55, stroke: 1pt, mark: (end: ">"))
    content((-1.1, 0.7), text(size: 9pt)[$u(t)$])
  }),
  caption: [Εναέριο όχημα σε περιστροφή γύρω από έναν άξονα. $phi(t)$ η γωνία προσανατολισμού, $u(t)$ η ροπή ελέγχου.],
) <fig-vehicle>

== Διατύπωση ως Γραμμικό Παραμετρικό Μοντέλο

Για την εφαρμογή της Μεθόδου Κλίσης, η @eq-plant γράφεται σε γραμμική παραμετρική μορφή $y(t) = theta^(*T) phi.alt(t)$, όπου ως έξοδος ορίζεται το μετρήσιμο σήμα $J dot.double(phi)(t)$ και ως διάνυσμα παλινδρόμησης το $phi.alt(t) = [dot(phi)(t), u(t)]^T$:

$ y(t) = J dot.double(phi)(t), quad
    phi.alt(t) = vec(dot(phi)(t), u(t)), quad
    theta^* = vec(-k, b) $ <eq-lpm>

Η επιλογή αυτή (για $d(t) = 0$) είναι δυνατή επειδή όλα τα σήματα που εμπλέκονται (συμπεριλαμβανομένου του $dot.double(phi)(t)$) θεωρούνται μετρήσιμα, οπότε δεν απαιτείται φιλτράρισμα των σημάτων με ευσταθές φίλτρο για την αποφυγή παραγώγων.

== Σύστημα Αναγνώρισης και Σφάλμα Εκτίμησης

Επιλέχθηκε σύστημα αναγνώρισης σε Παράλληλη Τοπολογία, δηλαδή χρησιμοποιώντας απευθείας τα μετρήσιμα σήματα του πραγματικού συστήματος ($dot(phi), u$) στον υπολογισμό της εκτιμώμενης εξόδου:

$ hat(y)(t) = hat(theta)^T (t) phi.alt(t) = -hat(k)(t) dot(phi)(t) + hat(b)(t) u(t) $ <eq-est>

όπου $hat(theta)(t) = [-hat(k)(t), hat(b)(t)]^T$ είναι η τρέχουσα εκτίμηση του $theta^*$. Το σφάλμα αναγνώρισης ορίζεται ως:

$ e(t) = y(t) - hat(y)(t)
    = J dot.double(phi)(t) + hat(k)(t) dot(phi)(t) - hat(b)(t) u(t) $ <eq-err>

και, ισοδύναμα με χρήση του παραμετρικού σφάλματος $tilde(theta) = hat(theta) - theta^*$, γράφεται $e = -tilde(theta)^T phi.alt$.

Η δομή αναγνώρισης απεικονίζεται στο @fig-identifier.

#figure(
  canvas(length: 1cm, {
    import cdraw: *

    // Plant block
    rect((0, 0.6), (3.2, 2.0), stroke: 1pt)
    content((1.6, 1.3), text(size: 9pt)[Σύστημα Αναγνώρισης])
    // arrow indicating parameter adaptation
    line((2.5, 2.05), (3.1, 2.65), mark: (end: ">"), stroke: 0.9pt)
    line((0.7, 0.25), (1, 0.55), stroke: 0.9pt)

    // Identifier block
    rect((0, -2.0), (3.2, -0.6), stroke: 1pt)
    content((1.6, -1.3), text(size: 9pt)[Πραγματικό Σύστημα])

    // Input u(t)
    line((-2.2, -1.3), (0, -1.3), mark: (end: ">"), stroke: 1pt)
    content((-1.85, -1), text(size: 9pt)[$u(t)$])
    // branch to identifier
    line((-1.1, 1.3), (-1.1, -1.3), stroke: 1pt)
    line((-1.1, 1.3), (0, 1.3), mark: (end: ">"), stroke: 1pt)

    // Output y = J phi-ddot
    line((3.2, 1.3), (5.2, 1.3), mark: (end: ">"), stroke: 1pt)
    content((5.6, 1.3), text(size: 9pt)[$hat(y)(t)$])

    // Estimated output
    line((3.2, -1.3), (5.2, -1.3), mark: (end: ">"), stroke: 1pt)
    content((5.6, -1.3), text(size: 9pt)[$y(t)$])

  }),
  caption: [Παράλληλη τοπολογία συστήματος αναγνώρισης για την εκτίμηση των $k$, $b$.],
) <fig-identifier>

Στην παρούσα εργασία επιλέχθηκε η παράλληλη δομή επειδή τα $phi, dot(phi), dot.double(phi)$ είναι όλα μετρήσιμα.

== Νόμος Προσαρμογής (Μέθοδος Κλίσης)

Η Μέθοδος Κλίσης αναζητά την ενημέρωση του $hat(theta)(t)$ που ελαχιστοποιεί μια τετραγωνική συνάρτηση κόστους του σφάλματος:

$ K(hat(theta)) = e^2 / 2 = (y - hat(theta)^T phi.alt)^2 / 2 $

Ο αναδρομικός νόμος προσαρμογής προκύπτει από την αρνητική κλίση της $K$ ως προς $hat(theta)$:

$ dot(hat(theta))(t) = -Gamma nabla K(hat(theta)) = Gamma e(t) phi.alt(t) $

όπου $Gamma = "diag"(gamma_1, gamma_2) > 0$ ο διαγώνιος πίνακας κερδών προσαρμογής. Αναλυτικά, για τις δύο συνιστώσες $hat(theta) = [-hat(k), hat(b)]^T$ και $phi.alt = [dot(phi), u]^T$, οι νόμοι ενημέρωσης γράφονται:

$ dot(hat(k))(t) & = -gamma_1 e(t) dot(phi)(t) \
  dot(hat(b))(t) & = +gamma_2 e(t) u(t) $ <eq-update>

Τα κέρδη $gamma_1, gamma_2 > 0$ ρυθμίζουν την ταχύτητα σύγκλισης. Μεγαλύτερες τιμές επιταχύνουν τη σύγκλιση, αλλά υπερβολικά μεγάλες καθιστούν το σύστημα δύσκολα επιλύσιμο αριθμητικά.

== Ανάλυση Ευστάθειας — Συνθήκη Επιμένουσας Διέγερσης

Θεωρώντας τη συνάρτηση Lyapunov $V = tilde(theta)^T Gamma^(-1) tilde(theta) slash 2 >= 0$ και παραγωγίζοντας ως προς τον χρόνο κατά μήκος των λύσεων της @eq-update προκύπτει:

$ dot(V) = tilde(theta)^T Gamma^(-1) dot(tilde(theta))
    = tilde(theta)^T e phi.alt
    = -e^2 <= 0 $

Συνεπώς $tilde(theta) in L_oo$ (τα παραμετρικά σφάλματα είναι φραγμένα) και $e in L_2$. Με εφαρμογή του Λήμματος του Barbalat προκύπτει $lim_(t -> oo) e(t) = 0$. Η ασυμπτωτική σύγκλιση των εκτιμήσεων στις πραγματικές τιμές, δηλαδή $hat(k)(t) -> k$ και $hat(b)(t) -> b$, απαιτεί επιπλέον να ικανοποιείται η συνθήκη επιμένουσας διέγερσης (ΣΕΔ) από το διάνυσμα παλινδρόμησης $phi.alt(t)$: να υπάρχουν σταθερές $alpha_0, alpha_1, T_0 > 0$ τέτοιες ώστε

$ alpha_1 I >= 1 / T_0 integral_t^(t+T_0) phi.alt(tau) phi.alt^T (tau) d tau >= alpha_0 I, quad forall t >= 0 $

Για την επιλεγμένη ημιτονοειδή είσοδο $u(t) = A sin(omega t)$ με $A = 0.25$, $omega = 0.5 pi$ rad/s, το σήμα $dot(phi)(t)$ στη μόνιμη κατάσταση ταλαντώνεται στην ίδια συχνότητα, καθώς η συνάρτηση μεταφοράς από $u$ στο $dot(phi)$ είναι $G(s) = b slash (J s + k)$. Συνεπώς $dot(phi)(t) = B sin(omega t + delta)$ με

$ B = (A b) / sqrt(k^2 + omega^2 J^2), quad delta = -arctan((omega J) / k) $

Υπολογίζοντας τον μέσο όρο του $phi.alt phi.alt^T$ σε μία περίοδο $T_0 = 2 pi slash omega$, προκύπτει:

$ M = 1 / T_0 integral_t^(t+T_0) phi.alt(tau) phi.alt^T (tau) d tau
    = 1 / 2 mat(B^2, A B cos delta; A B cos delta, A^2) $

Η ορίζουσα του $M$ είναι $det M = (A^2 B^2 slash 4) sin^2 delta$, η οποία είναι γνησίως θετική όσο $sin delta eq.not 0$. Και οι δύο ιδιοτιμές του $M$ είναι συνεπώς γνησίως θετικές, οπότε υπάρχουν $alpha_0, alpha_1 > 0$ με $alpha_0 I <= M <= alpha_1 I$ και η ΣΕΔ ικανοποιείται.

Αντικαθιστώντας τις αριθμητικές τιμές ($J = 0.025$, $k = 0.2$, $b = 1$), λαμβάνουμε $delta approx -0.194$ rad, $B approx 1.227$, και

$ M approx mat(0.7523, 0.1505; 0.1505, 0.0313), quad alpha_0 approx 1.1 dot 10^(-3), quad alpha_1 approx 0.783 $

Η μικρή τιμή του $alpha_0$ αναμένεται να οδηγήσει σε περιορισμένο ρυθμό σύγκλισης.

== Επίδραση Εξωτερικής Διαταραχής

Παρουσία διαταραχής $d(t) != 0$, η @eq-lpm γράφεται $y(t) = theta^(*T) phi.alt(t) + d(t)$ και το σφάλμα γίνεται $e = -tilde(theta)^T phi.alt + d$. Η παράγωγος της συνάρτησης Lyapunov προκύπτει:

$ dot(V) = -e^2 + e d $

η οποία δεν είναι πλέον αρνητικά ημιορισμένη, αλλά φραγμένη από $dot(V) <= -e^2 slash 2 + d^2 slash 2$. Το σύστημα εκτίμησης παραμένει ευσταθές, αλλά οι εκτιμήσεις $hat(k)(t), hat(b)(t)$ αναμένονται να μην συγκλίνουν πλέον στις πραγματικές τιμές, αλλά να εμφανίσουν μεροληψία και ταλαντώσεις.

== (α) Αποτελέσματα — Χωρίς Διαταραχή

Η προσομοίωση πραγματοποιήθηκε στο MATLAB με χρήση της `ode45`, για τις τιμές $k = 0.2$, $b = 1.0$, $J = 0.025$, κέρδη προσαρμογής $gamma_1 = gamma_2 = 200$, χρονικό ορίζοντα $T = 60$ s και αρχικές εκτιμήσεις $hat(k)(0) = 1.5 k = 0.3$, $hat(b)(0) = 1.5 b = 1.5$ (ώστε να είναι εμφανής η μεταβατική φάση σύγκλισης).

#figure(
  image("../outputs/phi_a.svg", width: 95%),
  caption: [Σύγκριση πραγματικής γωνίας $phi(t)$ και εκτιμώμενης $hat(phi)(t)$ από το παράλληλο μοντέλο, για $d(t) = 0$.],
) <fig-p1a-phi>

#figure(
  image("../outputs/err_a.svg", width: 95%),
  caption: [Εξέλιξη των σφαλμάτων $e_phi(t), e_k(t), e_b(t)$ χωρίς διαταραχή.],
) <fig-p1a-err>

Παρατηρείται η αναμενόμενη σύγκλιση των $e_k$, $e_b$ στο μηδέν σε περίπου 60 s. Το σφάλμα $e_phi$ του παράλληλου μοντέλου είναι της τάξης του $10^(-2)$ στη μόνιμη κατάσταση, πρακτικά αμελληταίο.

== (β) Αποτελέσματα — Με Εξωτερική Διαταραχή

Η ίδια προσομοίωση επαναλήφθηκε με $d(t) = 0.02 sin(2 t)$.

#figure(
  image("../outputs/phi_b.svg", width: 95%),
  caption: [Σύγκριση $phi(t)$ και $hat(phi)(t)$ παρουσία διαταραχής $d(t) = 0.02 sin(2t)$.],
) <fig-p1b-phi>

#figure(
  image("../outputs/err_b.svg", width: 95%),
  caption: [Εξέλιξη των σφαλμάτων $e_phi(t), e_k(t), e_b(t)$ με διαταραχή.],
) <fig-p1b-err>

Σε αντίθεση με την περίπτωση (α), τα σφάλματα δεν μηδενίζονται, αλλά εμφανίζουν ταλαντώσεις. Το σύστημα παρ'όλα αυτά 'φαίνεται να παραμένει φραγμένο.

= Θέμα 2

== Περιγραφή του Συστήματος

Θεωρούμε το δυναμικό σύστημα ενός απλού εκκρεμούς μήκους $l$ και μάζας $m$, το οποίο περιγράφεται από τη μη γραμμική εξίσωση:

$ dot.double(theta)(t) = -g / l sin(theta(t)) - c dot(theta)(t) + u(t) $ <eq-pend>

όπου $theta(t)$ [rad] είναι η γωνία του εκκρεμούς από την κατακόρυφο, $dot(theta)(t), dot.double(theta)(t)$ η γωνιακή ταχύτητα και επιτάχυνση, και $u(t)$ η εξωτερική ροπή ελέγχου. Η επιτάχυνση της βαρύτητας $g = 9.81 med "m/s"^2$ θεωρείται γνωστή, ενώ το μήκος $l > 0$ και ο συντελεστής απόσβεσης $c > 0$ είναι οι άγνωστες παράμετροι προς εκτίμηση. Για τα πειράματα επιλέγεται $u(t) = 0.5 sin(t)$ και υποθέτουμε ότι τα σήματα $theta(t), dot(theta)(t)$ και $u(t)$ είναι μετρήσιμα. Σε αντίθεση με το Θέμα 1, *το $dot.double(theta)(t)$ δεν θεωρείται μετρήσιμο*.

#figure(
  canvas(length: 1cm, {
    import cdraw: *

    // Ceiling
    line((-1.8, 0), (1.8, 0), stroke: 1pt)
    for i in range(9) {
      line((-1.6 + 0.4 * i, 0), (-1.85 + 0.4 * i, 0.28), stroke: 0.5pt)
    }

    // Pivot
    circle((0, 0), radius: 0.09, fill: black)

    // Vertical reference (dashed)
    line((0, 0), (0, -3.4), stroke: (dash: "dashed", thickness: 0.5pt, paint: luma(130)))

    // Rod at angle theta = 25 deg from vertical (leaning right)
    // x = L sin, y = -L cos, with L = 3
    let L = 3.0
    let s = 0.423 // sin 25
    let c_ = 0.906 // cos 25
    let bx = L * s
    let by = -L * c_
    line((0, 0), (bx, by), stroke: 1.4pt)

    // Bob
    circle((bx, by), radius: 0.32, fill: luma(220), stroke: 1.2pt)
    content((bx + 0.55, by - 0.05), text(size: 9pt)[$m$])

    // Length label
    content((bx / 2 + 0.35, by / 2), text(size: 9pt)[$l$])

    // Angle arc from vertical to rod
    arc((0, -0.95), start: -90deg, stop: -65deg, radius: 1.0, stroke: 0.9pt)
    content((0.28, -1.25), text(size: 10pt)[$theta$])

    // // Torque u(t) around pivot
    // arc((0, 0), start: 30deg, stop: 150deg, radius: 0.55, stroke: 1pt, mark: (end: ">"))
    // content((0, 0.92), text(size: 9pt)[$u(t)$])

    // Gravity arrow
    line((bx + 0.55, by + 1.75), (bx + 0.55, by + 1.05), mark: (end: ">"), stroke: 0.8pt, paint: luma(110))
    content((bx + 0.8, by + 1.55), text(size: 8pt, fill: luma(110))[$g$])
  }),
  caption: [Απλό εκκρεμές μήκους $l$ και μάζας $m$.],
) <fig-pendulum>

== Αναπαραμετροποίηση

Η @eq-pend είναι γραμμική ως προς τις σταθερές $g / l$ και $c$, αλλά μη γραμμική ως προς σκέτο $l$. Είναι επομένως βολικό να ορίσουμε τις αναπαραμετροποιημένες άγνωστες παραμέτρους:

$ theta_1^* = g / l, quad theta_2^* = c $

οπότε η δυναμική γράφεται:

$ dot.double(theta)(t) = -theta_1^* sin(theta(t)) - theta_2^* dot(theta)(t) + u(t) $

Ορίζοντας τις μεταβλητές κατάστασης $x_1 = theta, med x_2 = dot(theta)$, παίρνουμε το μη γραμμικό σύστημα:

$ dot(x)_1 & = x_2                                     \
  dot(x)_2 & = -theta_1^* sin(x_1) - theta_2^* x_2 + u $ <eq-pend-ss>

Μόλις εκτιμηθούν τα $hat(theta)_1, hat(theta)_2$, οι αρχικές παράμετροι μπορούν να ανακτηθούν ως $hat(l) = g / hat(theta)_1$ και $hat(c) = hat(theta)_2$.

== Σύστημα Αναγνώρισης

Σε αντίθεση με το Θέμα 1, όπου η διαθεσιμότητα του $dot.double(phi)$ επέτρεψε τη διατύπωση ενός αλγεβρικού γραμμικού παραμετρικού μοντέλου, εδώ το $dot.double(theta)$ δεν είναι μετρήσιμο, οπότε η Μέθοδος Κλίσης με απευθείας κατασκευή σφάλματος δεν εφαρμόζεται. Aντ' αυτού θα χρησιμοποιηθεί σχεδίαση κατά Lyapunov.

Ορίζουμε έναν *προσαρμοστικό observer 2ης τάξης* με έγχυση του σφάλματος θέσης $e_1$ και στις δύο εξισώσεις, και δική του κατάσταση $(hat(theta), hat(x)_2)$:

$ dot(hat(theta))(t) & = hat(x)_2 (t) + lambda_1 e_1 (t)                                                            \
  dot(hat(x))_2 (t)  & = -hat(theta)_1 (t) sin(theta(t)) - hat(theta)_2 (t) dot(theta)(t) + u(t) + lambda_2 e_1 (t) $ <eq-pend-est>

όπου $e_1 = theta - hat(theta)$ το σφάλμα θέσης και $(lambda_1, lambda_2) > 0$ τα κέρδη έγχυσης, τα οποία τοποθετούν τους πόλους του observer ως ρίζες του χαρακτηριστικού πολυωνύμου $s^2 + lambda_1 s + lambda_2 = 0$. Ορίζουμε επιπλέον το σφάλμα ταχύτητας $e_2 = dot(theta) - hat(x)_2$ και τα παραμετρικά σφάλματα $tilde(theta)_i = hat(theta)_i - theta_i^*$. Αφαιρώντας την @eq-pend-est από την @eq-pend-ss παίρνουμε τη δυναμική σφάλματος:

$ dot(e)_1 & = e_2 - lambda_1 e_1                                                    \
  dot(e)_2 & = -tilde(theta)_1 sin(theta) - tilde(theta)_2 dot(theta) - lambda_2 e_1 $ <eq-pend-err>

== Σχεδίαση κατά Lyapunov — Νόμοι Προσαρμογής

Θεωρούμε τη *σύνθετη* συνάρτηση Lyapunov:

$ V(e_1, e_2, tilde(theta)_1, tilde(theta)_2) = lambda_2 / 2 e_1^2 + 1 / 2 e_2^2 + 1 / (2 gamma_1) tilde(theta)_1^2 + 1 / (2 gamma_2) tilde(theta)_2^2 $

με $gamma_1, gamma_2 > 0$ κέρδη προσαρμογής. Το βάρος $lambda_2 slash 2$ μπροστά από το $e_1^2$ (αντί για $1 slash 2$) είναι κρίσιμο ώστε να ακυρωθούν οι όροι σύζευξης ($e_1 e_2$) στην παράγωγο. Παραγωγίζοντας κατά μήκος της @eq-pend-err:

$ dot(V) = lambda_2 e_1 dot(e)_1 + e_2 dot(e)_2 + 1 / gamma_1 tilde(theta)_1 dot(hat(theta))_1 + 1 / gamma_2 tilde(theta)_2 dot(hat(theta))_2 \
= -lambda_1 lambda_2 e_1^2 + underbrace(lambda_2 e_1 e_2 - lambda_2 e_1 e_2, =0) - e_2 tilde(theta)_1 sin(theta) - e_2 tilde(theta)_2 dot(theta) + 1 / gamma_1 tilde(theta)_1 dot(hat(theta))_1 + 1 / gamma_2 tilde(theta)_2 dot(hat(theta))_2 $

Επιλέγοντας τους νόμους προσαρμογής ώστε να ακυρωθούν οι εναπομείναντες όροι σύζευξης:

$ dot(hat(theta))_1 (t) & = gamma_1 e_2 (t) sin(theta(t)) \
  dot(hat(theta))_2 (t) & = gamma_2 e_2 (t) dot(theta)(t) $ <eq-pend-update>

προκύπτει:

$ dot(V) = -lambda_1 lambda_2 e_1^2 <= 0 $

Συνεπώς $e_1, e_2, tilde(theta)_1, tilde(theta)_2 in L_oo$ και $e_1 in L_2$· με εφαρμογή του Λήμματος Barbalat προκύπτει $lim_(t->oo) e_1(t) = 0$. Η *ασυμπτωτική σύγκλιση των παραμέτρων* $hat(theta)_i -> theta_i^*$ απαιτεί επιπλέον η είσοδος να ικανοποιεί την *συνθήκη επιμένουσας διέγερσης* ως προς το διάνυσμα παλινδρόμησης $phi.alt(t) = [sin(theta(t)), dot(theta)(t)]^T$. Για μονοσυχνοτική είσοδο $u(t) = 0.5 sin(t)$ και μικρές γωνίες ($sin theta approx theta$), τα δύο στοιχεία του $phi.alt$ είναι σχεδόν γραμμικώς εξαρτημένα (ολισθαίνουν φασικά κατά $90°$), οπότε *η ΣΕΔ δεν ικανοποιείται πλήρως* και αναμένεται ότι τα $hat(l), hat(c)$ θα παραμένουν φραγμένα αλλά δεν θα συγκλίνουν ασυμπτωτικά στις πραγματικές τιμές.

Το δομικό διάγραμμα του συστήματος αναγνώρισης φαίνεται στο @fig-identifier-pend.

#figure(
  canvas(length: 1cm, {
    import cdraw: *

    // Plant block
    rect((0, 0.6), (3.2, 2.0), stroke: 1pt)
    content((1.6, 1.3), text(size: 9pt)[Εκκρεμές])

    // Identifier block
    rect((0, -2.0), (3.2, -0.6), stroke: 1pt)
    content((1.6, -1.3), text(size: 9pt)[Σύστημα Αναγνώρισης])
    // parameter adaptation arrow
    line((2.7, -0.55), (3.1, -0.15), mark: (end: ">"), stroke: 0.9pt)

    // Input u(t)
    line((-2.2, 1.3), (0, 1.3), mark: (end: ">"), stroke: 1pt)
    content((-2.45, 1.3), text(size: 9pt)[$u(t)$])
    // branch to identifier
    line((-1.1, 1.3), (-1.1, -1.3), stroke: 1pt)
    line((-1.1, -1.3), (0, -1.3), mark: (end: ">"), stroke: 1pt)

    // Plant outputs x1, x2
    line((3.2, 1.55), (5.2, 1.55), mark: (end: ">"), stroke: 1pt)
    content((5.85, 1.55), text(size: 9pt)[$x_1, x_2$])
    // feedback of x1, x2 into identifier
    line((4.6, 1.55), (4.6, -1.6), stroke: 1pt)
    line((4.6, -1.6), (3.2, -1.6), mark: (end: ">"), stroke: 1pt)

    // Estimated x2
    line((3.2, -1.0), (5.2, -1.0), mark: (end: ">"), stroke: 1pt)
    content((5.75, -1.0), text(size: 9pt)[$hat(x)_2$])

    // Summing junction for e = x2 - x̂2
    circle((6.5, 0.3), radius: 0.25, stroke: 1pt)
    content((6.5, 0.3), text(size: 9pt)[$-$])
    // x2 into summer (from top branch)
    line((5.2, 1.55), (6.5, 1.55), stroke: 1pt)
    line((6.5, 1.55), (6.5, 0.55), mark: (end: ">"), stroke: 1pt)
    // x̂2 into summer
    line((5.2, -1.0), (6.5, -1.0), stroke: 1pt)
    line((6.5, -1.0), (6.5, 0.05), mark: (end: ">"), stroke: 1pt)

    // Error output
    line((6.75, 0.3), (7.6, 0.3), mark: (end: ">"), stroke: 1pt)
    content((7.85, 0.3), text(size: 9pt)[$e$])

    // Feedback of e to identifier (down and left)
    line((7.6, 0.3), (7.6, -2.6), stroke: 1pt)
    line((7.6, -2.6), (-1.7, -2.6), stroke: 1pt)
    line((-1.7, -2.6), (-1.7, -1.7), mark: (end: ">"), stroke: 1pt)
  }),
  caption: [Μεικτή δομή συστήματος αναγνώρισης για την εκτίμηση των $l$, $c$ του εκκρεμούς.],
) <fig-identifier-pend>

== Επίδραση Μετρητικού Θορύβου

Στο υπο-ερώτημα (β) οι μετρήσεις $theta(t), dot(theta)(t), u(t)$ υποτίθενται επιβαρυμένες με ημιτονοειδή θόρυβο $eta(t) = eta_0 sin(20 pi t)$. Αντικαθιστώντας τα θορυβώδη σήματα $theta_m = theta + eta, dot(theta)_m = dot(theta) + eta, u_m = u + eta$ στην @eq-pend-est, η δυναμική σφάλματος @eq-pend-err τροποποιείται ώστε να εμφανίζονται επιπλέον όροι εξαρτώμενοι από το $eta(t)$. Η παράγωγος Lyapunov παύει να είναι αρνητικά ημιορισμένη και παίρνει τη μορφή:

$ dot(V) = -lambda_1 lambda_2 e_1^2 + "όροι εξαρτώμενοι από " eta $

Το σύστημα παραμένει *ευσταθές υπό την έννοια της φραγμένης εισόδου – φραγμένης κατάστασης* (BIBS), αλλά οι εκτιμήσεις $hat(l), hat(c)$ εμφανίζουν:

- *Μόνιμη μεροληψία* (bias) λόγω της σύζευξης του θορύβου με τους όρους $sin(tilde(x)_1), tilde(x)_2$ στον νόμο προσαρμογής.
- *Ταλαντώσεις υψηλής συχνότητας* γύρω από τη μεροληπτική τιμή, με πλάτος που κλιμακώνεται με το $eta_0$ και τα κέρδη $gamma_1, gamma_2$.

Για την ποσοτικοποίηση της επίδρασης, θα παρουσιαστεί γραφική παράσταση του τελικού σφάλματος εκτίμησης $|e_l|, |e_c|$ συναρτήσει του πλάτους $eta_0$, ώστε να φανεί η (κατά προσέγγιση γραμμική) εξάρτηση της μεροληψίας από την ένταση του θορύβου.

== (α) Αποτελέσματα — Χωρίς Θόρυβο

Η προσομοίωση πραγματοποιήθηκε στο MATLAB με χρήση της `ode45`, για τις τιμές $l = 1.5$ m, $c = 0.5$, κέρδη προσαρμογής $gamma_1 = 5, gamma_2 = 2$, κέρδη έγχυσης observer $lambda_1 = 20, lambda_2 = 100$ (διπλός πόλος στο $s = -10$), χρονικό ορίζοντα $T = 40$ s, αρχική γωνία $theta(0) = 0.2$ rad και αρχικές εκτιμήσεις $hat(theta)_1 (0) = 0.9 dot g slash l approx 5.886$, $hat(theta)_2 (0) = 0.9 dot c = 0.45$ (δηλαδή $10%$ κάτω από τις αληθινές τιμές, προσομοιώνοντας ρεαλιστική αρχική αβεβαιότητα).

#figure(
  image("../outputs/theta_(a) no noise.svg", width: 95%),
  caption: [Σύγκριση πραγματικής γωνίας $theta(t)$ και εκτιμώμενης $hat(theta)(t)$ του observer, χωρίς θόρυβο.],
) <fig-p2a-theta>

#figure(
  image("../outputs/err_(a) no noise.svg", width: 95%),
  caption: [Εξέλιξη των σφαλμάτων $e_theta(t), e_l(t), e_c(t)$ χωρίς θόρυβο.],
) <fig-p2a-err>

*Σχόλια.* Ο observer παρακολουθεί το $theta(t)$ με υψηλή ακρίβεια (σφάλμα $|e_theta|_oo$ τάξης $10^(-3)$ rad $approx 0.06°$), επιβεβαιώνοντας την αποτελεσματικότητα της 2ης τάξης δομής έγχυσης. Τα σφάλματα παραμέτρων $e_l, e_c$ παραμένουν *φραγμένα* — όπως εγγυάται η ανάλυση Lyapunov — αλλά *δεν συγκλίνουν* στο μηδέν, παρουσιάζοντας αντίθετα μια αργή μονότονη ολίσθηση (π.χ. $e_l (T) approx +0.21$ m στα $T = 40$ s). Αυτή η συμπεριφορά είναι *θεωρητικά αναμενόμενη*: η μονοσυχνοτική είσοδος $u(t) = 0.5 sin(t)$ δεν εγείρει τη ΣΕΔ, οπότε ο νόμος προσαρμογής έχει μη-τετριμμένο μηδενόχωρο και δεν διακρίνει μεταξύ αλλαγών των $hat(theta)_1, hat(theta)_2$ που δίνουν ίδια έξοδο.

== (β) Αποτελέσματα — Με Μετρητικό Θόρυβο

Η προσομοίωση επαναλήφθηκε με θόρυβο $eta(t) = eta_0 sin(20 pi t)$ στις μετρήσεις $theta_m, dot(theta)_m, u_m$ για διάφορες τιμές πλάτους $eta_0 in {0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2}$, μετρώντας το τελικό σφάλμα εκτίμησης ως τον μέσο όρο του τελευταίου $1$ s της προσομοίωσης ώστε να έχουν προλάβει να συγκλίνουν οι εκτιμήσεις.

#figure(
  image("../outputs/err_vs_eta.svg", width: 95%),
  caption: [Τελικό σφάλμα εκτίμησης των παραμέτρων $|e_l|, |e_c|$ συναρτήσει του πλάτους θορύβου $eta_0$.],
) <fig-p2b-sweep>

Ο θόρυβος φαίνεται να επιδεινώνει ελαφρά την ακρίβεια για μικρά πλάτη θορύβου, αλλά δεν αλλάζει ποιοτικά τη συμπεριφορά, που κυριαρχείται από την αποτυχία της ΣΕΔ.

== Σύνοψη Θέματος 2

Η Μέθοδος Lyapunov παρέχει ισχυρές εγγυήσεις ευστάθειας. Τα σήματα παραμένουν φραγμένα ($cal(L)_oo$), το σφάλμα κατάστασης $e_1 -> 0$ και το σύστημα είναι ανθεκτικό σε υψίσυχνο θόρυβο. Ωστόσο, χωρίς να ισχύει η ΣΕΔ δεν είανι σίγουρη η σύγκλιση των παραμέτρων. Με τη ημιτονοειδή είσοδο $u(t) = 0.5 sin(t)$ τα $e_l, e_c$ παρουσιάζουν αργή ολίσθηση.
