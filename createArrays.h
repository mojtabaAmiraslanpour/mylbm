    // Defining rho
    double** rho_ = new double* [lx_ + 2];
    for (int i = 0; i < lx_ + 2; i++) {
        rho_[i] = new double[ly_ + 2];
        for (int j = 0; j < ly_ + 2; j++) {
            rho_[i][j] = 1.0;
        }
    }

    // Defining Ux
    double y_;
    int L = ly_ - 2;
    double** Ux_ = new double* [lx_ + 2];
    for (int i = 0; i < lx_ + 2; i++) {
        Ux_[i] = new double[ly_ + 2];
        for (int j = 0; j < ly_ + 2; j++) {
            y_ = j - 1.5;
            Ux_[i][j] = 4. * uMax_ / (L * L) * (L * y_ - y_ * y_);
        }
    }

    // Defining Uy
    double** Uy_ = new double* [lx_ + 2];
    for (int i = 0; i < lx_ + 2; i++) {
        Uy_[i] = new double[ly_ + 2];
        for (int j = 0; j < ly_ + 2; j++) {
            Uy_[i][j] = 0;
        }
    }

    // Defining population array f(q,lx,ly)
    double*** f_;
    f_ = new double** [q_];
    for (int k = 0; k < q_; k++) {
        f_[k] = new double* [lx_ + 2];
        for (int i = 0; i < lx_ + 2; i++) {
            f_[k][i] = new double[ly_ + 2];
            for (int j = 0; j < ly_ + 2; j++) {
                f_[k][i][j] = 0;
            }
        }
    }

    // Defining population array fStar(q,lx,ly)
    double*** fStar_;
    fStar_ = new double** [q_];
    for (int k = 0; k < q_; k++) {
        fStar_[k] = new double* [lx_ + 2];
        for (int i = 0; i < lx_ + 2; i++) {
            fStar_[k][i] = new double[ly_ + 2];
            for (int j = 0; j < ly_ + 2; j++) {
                fStar_[k][i][j] = 0;
            }
        }
    }

    // Defining population array fEq(q,lx,ly)
    double*** fEq_;
    fEq_ = new double** [q_];
    for (int k = 0; k < q_; k++) {
        fEq_[k] = new double* [lx_ + 2];
        for (int i = 0; i < lx_ + 2; i++) {
            fEq_[k][i] = new double[ly_ + 2];
            for (int j = 0; j < ly_ + 2; j++) {
                fEq_[k][i][j] = 0;
            }
        }
    }

    // // Defining c * u
    // double*** cu_;
    // cu_ = new double** [q_];
    // for (int k = 0; k < q_; k++) {
    //     cu_[k] = new double* [lx_ + 2];
    //     for (int i = 0; i < lx_ + 2; i++) {
    //         cu_[k][i] = new double[ly_ + 2];
    //         for (int j = 0; j < ly_ + 2; j++) {
    //             cu_[k][i][j] = 0;
    //         }
    //     }
    // }

    // Defining c * u
    double** cu_ = new double* [lx_ + 2];
    for (int i = 0; i < lx_ + 2; i++) {
        cu_[i] = new double[ly_ + 2];
        for (int j = 0; j < ly_ + 2; j++) {
            cu_[i][j] = 0;
        }
    }